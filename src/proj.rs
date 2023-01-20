//!
//! Projection installation
//!

use crate::datum_params::DatumParams;
use crate::datum_transform::Datum;
use crate::datums::{self, DatumDefn};
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::projections::{ProjFn, ProjParams};
use crate::{ellipsoids, prime_meridians, projstring, units};

use std::fmt;

pub type Axis = [u8; 3];

const NORMALIZED_AXIS: Axis = [b'e', b'n', b'u'];

/// A Proj obect hold informations and parameters
/// for a projection
pub struct Proj {
    pub(crate) pm: f64,
    pub(crate) ellps: Ellipsoid,
    pub(crate) datum: Datum,
    pub(crate) axis: Axis,
    pub(crate) is_geocent: bool,
    pub(crate) is_latlong: bool,
    pub(crate) to_meter: f64,
    pub(crate) vto_meter: f64,
    pub(crate) x0: f64,
    pub(crate) y0: f64,
    pub(crate) k0: f64,
    pub(crate) lam0: f64,
    pub(crate) phi0: f64,
    pub(crate) geoc: bool,
    pub(crate) over: bool, // over-ranging flag
    pub(crate) projname: &'static str,
    pub(crate) projdata: ProjParams,
    pub(crate) inverse: Option<ProjFn>,
    pub(crate) forward: Option<ProjFn>,
}

impl Proj {
    // ----------------
    // Datum definition
    // ----------------
    fn datum_defn<'a>(params: &'a ParamList) -> Result<Option<&'a DatumDefn>> {
        // Do we have a "datum" parameter ?
        params
            .get("datum")
            .map(|p| match datums::find_datum(p.try_into()?) {
                Some(v) => Ok(Some(v)),
                None => Err(Error::InvalidDatum),
            })
            .unwrap_or(Ok(None))
    }

    // --------------
    // Prime meridian
    // --------------
    fn prime_meridian(params: &ParamList) -> Result<f64> {
        params
            .get("pm")
            .map(
                |p| match prime_meridians::find_prime_meridian(p.try_into()?) {
                    Some(v) => Ok(v),
                    None => p.try_into(),
                },
            )
            .unwrap_or(Ok(0.))
    }

    // -----------------
    // Datum parameters
    // ----------------
    fn datum_params(params: &ParamList, defn: Option<&DatumDefn>) -> Result<DatumParams> {
        // Precedence order is 'nadgrids', 'towgs84', 'datum'
        if let Some(p) = params.get("nadgrids") {
            // Nadgrids
            DatumParams::from_nagrid_str(p.try_into()?)
        } else if let Some(p) = params.get("towgs84") {
            DatumParams::from_towgs84_str(p.try_into()?)
            // ToWGS84
        } else if let Some(p) = defn {
            DatumParams::try_from(&p.params)
        } else {
            Ok(DatumParams::default())
        }
    }

    // -----------------
    // Ellipsoid
    // ----------------
    fn ellipsoid(params: &ParamList, datum_def: Option<&DatumDefn>) -> Result<Ellipsoid> {
        if let Some(radius) = params.get("R") {
            // Sphere override everything
            Ellipsoid::sphere(radius.try_into()?)
        } else if let Some(p) = params.get("ellps") {
            // Return from ellipse definition
            match ellipsoids::find_ellipsoid(p.try_into()?) {
                Some(defn) => Ellipsoid::try_from_ellipsoid_with_params(defn, params),
                None => Err(Error::InvalidEllipsoid),
            }
        } else if let Some(defn) = datum_def {
            // Retrieve from datum definition + parameters
            Ellipsoid::try_from_ellipsoid_with_params(defn.ellps, params)
        } else {
            // Get a free WGS84
            Ellipsoid::try_from_ellipsoid_with_params(&ellipsoids::constants::WGS84, params)
        }
    }

    // -----------------
    // Axis
    // ----------------
    fn axis(params: &ParamList) -> Result<Axis> {
        if let Some(p) = params.get("axis") {
            let axis_arg: &str = p.try_into()?;
            if axis_arg.len() != 3 {
                Err(Error::InvalidAxis)
            } else {
                let mut axis = [0u8, 0u8, 0u8];
                // Find Easting/Westing
                // This ensure that no token is repeated unless
                // one of the `find` will fail.
                let ew = axis_arg.find(['e', 'w']).ok_or(Error::InvalidAxis)?;
                let ns = axis_arg.find(['n', 's']).ok_or(Error::InvalidAxis)?;
                let ud = axis_arg.find(['u', 'd']).ok_or(Error::InvalidAxis)?;
                axis[ew] = axis_arg.as_bytes()[ew];
                axis[ns] = axis_arg.as_bytes()[ns];
                axis[ud] = axis_arg.as_bytes()[ud];
                Ok(axis)
            }
        } else {
            Ok(NORMALIZED_AXIS)
        }
    }

    /// Return true if the axis are normalized
    pub fn normalized_axis(&self) -> bool {
        self.axis == NORMALIZED_AXIS
    }

    // -----------------
    // Units
    // ----------------
    fn units(params: &ParamList, name: &str, default: f64) -> Result<f64> {
        if let Some(p) = params.get(name) {
            units::find_unit_to_meter(p.try_into()?)
                .map(Ok)
                .unwrap_or_else(|| p.try_into())
        } else {
            Ok(default)
        }
    }

    /// Consume a ParamList and create a Proj object
    pub fn init(params: ParamList) -> Result<Self> {
        // Projection name
        let _projname = params.get("proj").ok_or(Error::MissingProjectionError);

        // Get datum definition (if any)
        let datum_defn = Self::datum_defn(&params)?;

        // Get datum parameters
        let datum_params = Self::datum_params(&params, datum_defn)?;

        // Do we have an ellipse ?
        let ellps = Self::ellipsoid(&params, datum_defn)?;

        // Get prime meridian
        let pm = Self::prime_meridian(&params)?;

        // Axis
        let axis = Self::axis(&params)?;

        // units
        let to_meter = Self::units(&params, "to_meter", 1.)?;
        // XXX in proj4 vto_meter accept fractional expression: '/'
        let vto_meter = Self::units(&params, "vto_meter", to_meter)?;

        // Datum
        let datum = Datum::new(&ellps, datum_params);

        Ok(Self {
            pm,
            ellps,
            datum,
            axis,
            is_geocent: false,
            is_latlong: false,
            to_meter,
            vto_meter,
            // Central meridian_
            lam0: params.try_value("lon_0", 0.)?,
            phi0: params.try_value("lat_0", 0.)?,
            x0: params.try_value("x_0", 0.)?,
            y0: params.try_value("y_0", 0.)?,
            // Proj4 compatibility
            k0: params
                .get("k0")
                .or_else(|| params.get("k"))
                .map(|p| p.try_into())
                .unwrap_or(Ok(1.))?,
            geoc: false,
            over: params.check_option("over")?,
            projdata: ProjParams::NoParams,
            projname: "",
            inverse: None,
            forward: None,
        }
        .prepare())
    }

    fn prepare(self) -> Self {
        self
    }

    /// Create from projstring definition
    pub fn from_projstr(s: &str) -> Result<Self> {
        Self::init(projstring::parse(s)?)
    }

    /// Create projection from string
    pub fn from_user_string(s: &str) -> Result<Self> {
        let s = s.trim();
        if s.starts_with('+') {
            Self::from_projstr(s)
        } else if s.eq_ignore_ascii_case("WGS84") {
            Self::from_projstr("+proj=longlat +ellps=WGS84")
        } else {
            Err(Error::UnrecognizedFormat)
        }
    }
}

// -------------
// Display
// -------------
impl fmt::Debug for Proj {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "pm:         {:?}", self.pm)?;
        write!(f, "ellps:      {:#?}", self.ellps)?;
        write!(f, "datum:      {:#?}", self.datum)?;
        write!(f, "axis:       {:?}", self.axis)?;
        write!(f, "is_geocent: {:?}", self.is_geocent)?;
        write!(f, "is_latlong: {:?}", self.is_latlong)?;
        write!(f, "to_meter:   {:?}", self.to_meter)?;
        write!(f, "vto_meter:  {:?}", self.vto_meter)?;
        write!(f, "x0:         {:?}", self.x0)?;
        write!(f, "y0:         {:?}", self.y0)?;
        write!(f, "lam0        {:?}", self.lam0)?;
        write!(f, "phi0:       {:?}", self.phi0)?;
        write!(f, "geoc:       {:?}", self.geoc)?;
        write!(f, "over:       {:?}", self.over)?;
        write!(f, "projname:   {:?}", self.projname)?;
        write!(f, "projdata:   {:#?}", self.projdata)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use crate::errors::{Error, Result};

    const EPSG_102018: &str = concat!(
        "+proj=gnom +lat_0=90 +lon_0=0 +x_0=6300000 +y_0=6300000",
        "+ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    );

    const TESTMERC: &str = "+proj=merc +lon_0=5.937 +lat_ts=45.027 +ellps=sphere";
    const TESTMERC2: &str = concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 ",
        "+units=m +k=1.0 +nadgrids=@null +no_defs"
    );
    const INVALID_ELLPS: &str = "+proj=merc +lon_0=5.937 +lat_ts=45.027 +ellps=foo";

    #[test]
    fn proj_test_EPSG_102018() {
        let p: Proj = Proj::from_projstr(EPSG_102018).unwrap();
    }

    #[test]
    fn proj_test_merc() {
        let p: Proj = Proj::from_projstr(TESTMERC).unwrap();
    }

    #[test]
    fn proj_test_merc2() {
        let p: Proj = Proj::from_projstr(TESTMERC2).unwrap();
    }

    #[test]
    fn proj_invalid_ellps_param() {
        let p: Result<Proj> = Proj::from_projstr(INVALID_ELLPS);

        assert!(p.is_err());
        let err = p.unwrap_err();
        assert!(matches!(err, Error::InvalidEllipsoid));
    }
}
