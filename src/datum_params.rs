//!
//! Handle datum parameters
//!
use crate::constants::SEC_TO_RAD;
use crate::datums::DatumParamDefn;
use crate::errors::{Error, Result};
use crate::nadgrids::NadgridShift;

/// Datum parameters
#[derive(Clone, Debug, PartialEq)]
pub enum DatumParams<GS: NadgridShift> {
    ToWGS84_0,
    ToWGS84_3(f64, f64, f64),
    ToWGS84_7(f64, f64, f64, f64, f64, f64, f64),
    NadGrids(GS),
    NoDatum,
}

impl<GS: NadgridShift> Default for DatumParams<GS> {
    fn default() -> Self {
        DatumParams::NoDatum
    }
}

impl<GS: NadgridShift> DatumParams<GS> {
    /// Create parameters from a 'towgs84 like string'
    /// Values are expected to be in second of arcs
    pub fn from_towgs84_str(towgs84: &str) -> Result<Self> {
        let mut i = towgs84.split(',');

        #[inline]
        fn parse(v: Option<&str>) -> Result<f64> {
            v.unwrap()
                .trim()
                .parse::<f64>()
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

    pub fn from_nagrid_str(nadgrids: &str) -> Result<Self> {
        GS::new_grid_transform(nadgrids).map(|g| Self::NadGrids(g))
    }
}

// Convert from datum parameters definition
impl<GS: NadgridShift> TryFrom<&DatumParamDefn> for DatumParams<GS> {
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
            DatumParamDefn::NadGrids(s) => Self::from_nagrid_str(s),
        }
    }
}
