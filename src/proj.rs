//!
//! Projection installation
//!

use crate::datum_params::DatumParams;
use crate::datum_transform::Datum;
use crate::datums::{self, DatumDefn};
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::nadgrids::{NadgridShift, NullGridShift};
use crate::parameters::ParamList;
use crate::{ellipsoids, prime_meridians, projstring, units};
//============================
//
// Projection
//
//============================

pub type Axis = [u8; 3];

const NORMALIZED_AXIS: Axis = [b'e', b'n', b'u'];

#[derive(Debug)]
pub struct Projection<N: NadgridShift = NullGridShift> {
    pub(crate) pm: f64,
    pub(crate) ellps: Ellipsoid,
    pub(crate) datum: Datum<N>,
    pub(crate) axis: Axis,
    pub(crate) to_meter: f64,
    pub(crate) is_geocent: bool,
    pub(crate) is_latlong: bool,
}

impl<N: NadgridShift> Projection<N> {
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
                    None => p.try_convert::<f64>(),
                },
            )
            .unwrap_or(Ok(0.))
    }

    // -----------------
    // Datum parameters
    // ----------------
    fn datum_params(params: &ParamList, defn: Option<&DatumDefn>) -> Result<DatumParams<N>> {
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
        return self.axis == NORMALIZED_AXIS;
    }

    // -----------------
    // Units
    // ----------------
    fn units(params: &ParamList) -> Result<f64> {
        if let Some(p) = params.get("to_meter") {
            units::find_unit_to_meter(p.try_into()?)
                .map(Ok)
                .unwrap_or_else(|| p.try_convert::<f64>())
        } else {
            Ok(1.)
        }
    }

    /// Prepare the projection
    fn prepare(self) -> Self {
        self
    }

    /// Consume a ParamList and create a Projection object
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
        let to_meter = Self::units(&params)?;

        // Datum
        let datum = Datum::new(&ellps, datum_params);

        Ok(Self {
            pm,
            ellps,
            datum,
            axis,
            to_meter,
            is_geocent: false,
            is_latlong: false,
        }
        .prepare())
    }

    /// Create projection from string
    pub fn from_projstr(s: &str) -> Result<Self> {
        Self::init(projstring::parse(s)?)
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
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0",
        "+units=m +k=1.0 +nadgrids=@null +no_defs"
    );
    const INVALID_ELLPS: &str = "+proj=merc +lon_0=5.937 +lat_ts=45.027 +ellps=foo";

    #[test]
    fn proj_test_EPSG_102018() {
        let p: Projection = Projection::from_projstr(EPSG_102018).unwrap();
    }

    #[test]
    fn proj_test_merc() {
        let p: Projection = Projection::from_projstr(TESTMERC).unwrap();
    }

    #[test]
    fn proj_test_merc2() {
        let p: Projection = Projection::from_projstr(TESTMERC2).unwrap();
    }

    #[test]
    fn proj_invalid_ellps_param() {
        let p: Result<Projection> = Projection::from_projstr(INVALID_ELLPS);

        assert!(p.is_err());
        assert!(matches!(p.unwrap_err(), Error::InvalidEllipsoid));
    }
}
