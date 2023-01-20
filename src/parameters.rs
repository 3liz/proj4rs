//!
//! Projection parameters
//!
//!
use crate::errors::{Error, Result};
use std::str::FromStr;

// XXX Parsing code take about 1kb in wasm, try to use JS parsing functions
// parseInt(), parseFloat() with wasm-bingen


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

impl<'a> Parameter<'a> {
    fn try_value<F: FromStr>(&self) -> Result<F> {
        match self.value.map(F::from_str) {
            None => Err(Error::NoValueParameter(self.name.into())),
            Some(result) => result.map_err(|_err| Error::ParameterValueError(self.name.into())),
        }
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
    /// Create a Parameter list from a vector of Params
    pub fn new(params: Vec<Parameter<'a>>) -> Self {
        Self(params)
    }

    /// Return Some(param) if the parameter `name` exists `None` otherwise.
    pub fn get(&self, name: &str) -> Option<&Parameter<'a>> {
        self.0.iter().find(|p| p.name == name)
    }

    pub fn check_option(&self, name: &str) -> Result<bool> {
        self.get(name)
            .map(|p| p.check_option())
            .unwrap_or(Ok(false))
    }

    pub fn try_value<T>(&self, name: &str, default: T) -> Result<T>
    where
        T: FromStr,
    {
        self.get(name)
            .map(|p| p.try_value::<T>())
            .unwrap_or(Ok(default))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

        assert_eq!(params.try_value::<f64>("foo", 0.).unwrap(), 1234.);
        assert_eq!(params.try_value::<f64>("bar", 0.).unwrap(), 0.);
    }
}
