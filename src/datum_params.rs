//!
//! Handle datum parameters
//!
use crate::datums::DatumParamDefn;
use crate::errors::{Error, Result};
use crate::math::consts::SEC_TO_RAD;
use crate::nadgrids::NadGrids;
use crate::parse::FromStr;

/// Datum parameters
#[derive(Debug, PartialEq)]
pub enum DatumParams {
    ToWGS84_0,
    ToWGS84_3(f64, f64, f64),
    ToWGS84_7(f64, f64, f64, f64, f64, f64, f64),
    NadGrids(NadGrids),
    NoDatum,
}

impl Default for DatumParams {
    fn default() -> Self {
        DatumParams::NoDatum
    }
}

impl DatumParams {
    /// Create parameters from a 'towgs84 like string'
    /// Values are expected to be in second of arcs
    pub fn from_towgs84_str(towgs84: &str) -> Result<Self> {
        let mut i = towgs84.split(',');

        // XXX Use js_sys::parsefloat with Wasm
        // It save about 20ko !
        fn parse(v: Option<&str>) -> Result<f64> {
            f64::from_str(v.unwrap_or("").trim())
                .map(|v| v * SEC_TO_RAD)
                .map_err(|_| Error::InvalidToWGS84String)
        }

        match towgs84.split(',').count() {
            3 => Ok(DatumParams::ToWGS84_3(
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
            )),
            7 => Ok(DatumParams::ToWGS84_7(
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
                parse(i.next())?,
            )),
            _ => Err(Error::InvalidToWGS84String),
        }
    }

    pub fn from_nadgrid_str(nadgrids: &str) -> Result<Self> {
        NadGrids::new_grid_transform(nadgrids).map(Self::NadGrids)
    }

    pub fn use_nadgrids(&self) -> bool {
        matches!(self, Self::NadGrids(_))
    }

    pub fn no_datum(&self) -> bool {
        matches!(self, Self::NoDatum)
    }
}

// Convert from datum parameters definition
impl TryFrom<&DatumParamDefn> for DatumParams {
    type Error = Error;

    fn try_from(defn: &DatumParamDefn) -> Result<Self> {
        match defn {
            DatumParamDefn::ToWGS84_0 => Ok(Self::ToWGS84_0),
            DatumParamDefn::ToWGS84_3(dx, dy, dz) => Ok(Self::ToWGS84_3(*dx, *dy, *dz)),
            DatumParamDefn::ToWGS84_7(dx, dy, dz, rx, ry, rz, s) => Ok(Self::ToWGS84_7(
                *dx,
                *dy,
                *dz,
                *rx * SEC_TO_RAD,
                *ry * SEC_TO_RAD,
                *rz * SEC_TO_RAD,
                *s / 1_000_000.0 + 1.,
            )),
            DatumParamDefn::NadGrids(s) => Self::from_nadgrid_str(s),
        }
    }
}
