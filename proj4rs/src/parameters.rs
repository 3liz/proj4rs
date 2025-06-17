//!
//! Projection parameters
//!
//!
use crate::errors::{Error, Result};
use crate::parse::FromStr;

/// Struct holding a pair key/value
pub struct Parameter<'a> {
    pub name: &'a str,
    pub value: Option<&'a str>,
}

impl<'a> TryFrom<&Parameter<'a>> for f64 {
    type Error = Error;

    fn try_from(p: &Parameter<'a>) -> Result<f64> {
        p.try_value::<f64>()
    }
}

// Param -> i32
impl<'a> TryFrom<&Parameter<'a>> for i32 {
    type Error = Error;

    fn try_from(p: &Parameter<'a>) -> Result<i32> {
        p.try_value::<i32>()
    }
}

// Param -> &str
impl<'a> TryFrom<&Parameter<'a>> for &'a str {
    type Error = Error;

    fn try_from(p: &Parameter<'a>) -> Result<&'a str> {
        p.value
            .ok_or_else(|| Error::NoValueParameter(p.name.into()))
    }
}

impl Parameter<'_> {
    fn try_value<F: FromStr>(&self) -> Result<F> {
        match self.value.map(F::from_str) {
            None => Err(Error::NoValueParameter(self.name.into())),
            Some(result) => result.map_err(|_| Error::ParameterValueError(self.name.into())),
        }
    }

    /// Return a value in radians
    ///
    /// By default it is assumed that unspecified input
    /// is in degree.
    ///
    /// If value string end with 'r', then it assumed that
    /// input value is given directly in radians
    ///
    /// Syntax suffix is [rR][EeWwNnSs]
    fn try_angular_value(&self) -> Result<f64> {
        const FWD_SFX: &[char; 4] = &['E', 'e', 'N', 'n'];
        const INV_SFX: &[char; 4] = &['W', 'w', 'S', 's'];
        const RAD_SFX: &[char; 2] = &['r', 'R'];

        match self.value {
            None => Err(Error::NoValueParameter(self.name.into())),
            Some(s) => {
                let (s, sgn) = if let Some(s) = s.strip_suffix(FWD_SFX) {
                    (s, 1.0)
                } else if let Some(s) = s.strip_suffix(INV_SFX) {
                    (s, -1.0)
                } else {
                    (s, 1.0)
                };

                if let Some(s) = s.strip_suffix(RAD_SFX) {
                    f64::from_str(s).map(|v| sgn * v)
                } else {
                    // Degrees
                    Self::parse_dms(s).map(|v| sgn * v.to_radians())
                }
            }
            .map_err(|_| Error::ParameterValueError(self.name.into())),
        }
    }

    /// DMS input parser
    ///
    /// DMS input is expected to be DD[dD]MM'SS"
    fn parse_dms(s: &str) -> Result<f64, <f64 as FromStr>::Err> {
        const DEG_SFX: &[char; 3] = &['d', 'D', '\u{00b0}'];

        fn parse_number_part<'a>(
            s: &'a str,
            sfx: &'static str,
        ) -> Result<(&'a str, f64), <f64 as FromStr>::Err> {
            // check the first non-numeric symbol starting from
            // the end and return prefix and parsed value.
            Ok(if let Some(s) = s.strip_suffix(sfx) {
                if let Some((pfx, v)) = s.rsplit_once(|c: char| !c.is_ascii_digit() && c != '.') {
                    (&s[..(pfx.len() + 1)], f64::from_str(v)?)
                } else {
                    ("", f64::from_str(s)?)
                }
            } else {
                (s, 0.)
            })
        }

        if s.is_empty() {
            // Force error if string is empty
            return f64::from_str(s);
        }

        let (s, seconds) = parse_number_part(s, "\"")?;
        let (s, minutes) = parse_number_part(s, "'")?;

        let s = s.trim_end_matches(DEG_SFX);
        let degrees = if !s.is_empty() { f64::from_str(s)? } else { 0. };

        Ok(degrees + (minutes + seconds / 60.) / 60.)
    }

    /// Check the token as a boolean flag
    ///
    /// Return true if the token is present alone (no value), false
    /// if the token is not present or parse the
    /// value as bool if any (either 'true' or 'false')
    pub fn check_option(&self) -> Result<bool> {
        self.value
            .map(bool::from_str)
            .unwrap_or(Ok(true))
            .map_err(|_err| Error::ParameterValueError(self.name.into()))
    }
}

/// List of parameters
pub struct ParamList<'a>(Vec<Parameter<'a>>);

impl<'a> ParamList<'a> {
    /// Return Some(param) if the parameter `name` exists `None` otherwise.
    pub fn get(&self, name: &str) -> Option<&Parameter<'a>> {
        self.0.iter().find(|p| p.name == name)
    }

    pub fn check_option(&self, name: &str) -> Result<bool> {
        self.get(name)
            .map(|p| p.check_option())
            .unwrap_or(Ok(false))
    }

    pub fn try_value<T>(&'a self, name: &str) -> Result<Option<T>>
    where
        T: TryFrom<&'a Parameter<'a>, Error = Error>,
    {
        self.get(name).map(|p| T::try_from(p)).transpose()
    }

    pub fn try_angular_value(&self, name: &str) -> Result<Option<f64>> {
        self.get(name).map(|p| p.try_angular_value()).transpose()
    }
}

// Create from Parameter iterator
impl<'a> FromIterator<Parameter<'a>> for ParamList<'a> {
    fn from_iter<I: IntoIterator<Item = Parameter<'a>>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

#[cfg(test)]
mod tests {
    use super::Parameter;
    use crate::projstring::parse;

    #[test]
    fn param_check_options() {
        let params = parse("+foo +bar=true +baz=false +bad=foobar").unwrap();

        assert_eq!(params.check_option("foo").unwrap(), true);
        assert_eq!(params.check_option("bar").unwrap(), true);
        assert_eq!(params.check_option("baz").unwrap(), false);
        assert_eq!(params.check_option("foobar").unwrap(), false);

        assert!(params.check_option("bad").is_err());
    }

    #[test]
    fn param_int_as_float() {
        let params = parse("+foo=0 +bar=1234 +baz=-2").unwrap();

        assert_eq!(f64::try_from(params.get("foo").unwrap()).unwrap(), 0.);
        assert_eq!(f64::try_from(params.get("bar").unwrap()).unwrap(), 1234.);
        assert_eq!(f64::try_from(params.get("baz").unwrap()).unwrap(), -2.);
    }

    #[test]
    fn param_try_value() {
        let params = parse("+foo=1234").unwrap();

        assert_eq!(params.try_value::<f64>("foo").unwrap().unwrap_or(0.), 1234.);
        assert_eq!(params.try_value::<f64>("bar").unwrap().unwrap_or(0.), 0.);
    }

    #[test]
    fn param_try_angular_value() {
        // Invalid specifier
        let params = parse("+foo").unwrap();
        assert!(params.try_angular_value("foo").is_err());

        // Input in degrees
        let params = parse("+foo=47").unwrap();
        assert_eq!(
            params.try_angular_value("foo").unwrap().unwrap(),
            47.0f64.to_radians()
        );

        // Input in radians
        let params = parse("+foo=2.3r").unwrap();
        assert_eq!(params.try_angular_value("foo").unwrap().unwrap(), 2.3f64);

        let params = parse("+foo=2.3R").unwrap();
        assert_eq!(params.try_angular_value("foo").unwrap().unwrap(), 2.3f64);

        // Orientation specifier
        let params = parse("+foo=47w").unwrap();
        assert_eq!(
            params.try_angular_value("foo").unwrap().unwrap(),
            -47.0f64.to_radians()
        );

        let params = parse("+foo=47W").unwrap();
        assert_eq!(
            params.try_angular_value("foo").unwrap().unwrap(),
            -47.0f64.to_radians()
        );

        // Input in radians with orientation specifier
        let params = parse("+foo=2.3rw").unwrap();
        assert_eq!(params.try_angular_value("foo").unwrap().unwrap(), -2.3);

        let params = parse("+foo=2.3R").unwrap();
        assert_eq!(params.try_angular_value("foo").unwrap().unwrap(), 2.3);

        // Invalid specifier
        let params = parse("+foo=2.3wr").unwrap();
        assert!(params.try_angular_value("foo").is_err());

        // DWS value
        let params = parse("+foo=38d30'9\"").unwrap();
        assert_eq!(
            params.try_angular_value("foo").unwrap().unwrap(),
            38.5025_f64.to_radians(),
        );
    }

    #[test]
    fn param_dms_parsing() {
        assert_eq!(Parameter::parse_dms("38"), Ok(38.0));
        assert_eq!(Parameter::parse_dms("38d"), Ok(38.0));
        assert_eq!(Parameter::parse_dms("38D"), Ok(38.0));
        assert_eq!(Parameter::parse_dms("38\u{00b0}"), Ok(38.0));
        assert_eq!(Parameter::parse_dms("38d30'"), Ok(38.5));
        assert_eq!(Parameter::parse_dms("38d30'9\""), Ok(38.5025));
        assert_eq!(Parameter::parse_dms("38d30.15'"), Ok(38.5025));
        assert_eq!(Parameter::parse_dms("30'9\""), Ok(0.5025));
    }
}
