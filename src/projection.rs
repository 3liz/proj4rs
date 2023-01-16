//!
//! Derived values from projection definition.
//!
//! Precomputed for optimization
//!

use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::nadgrids::{NadgridProvider, NadgridShift};
use crate::parameters::ParamList;
use crate::{datums, ellipsoids, prime_meridians, projstring};

pub struct Projection {
    pj_pm: f64,
    pj_ellps: Ellipsoid,
}

impl Projection {
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

    /// Consume a  ParamList and create a Projection object
    pub fn try_new(params: ParamList) -> Result<Self> {
        // Projection name
        let projname = params.get("proj").ok_or(Error::MissingProjectionError);

        // Get prime meridian
        let pj_pm = Self::prime_meridian(&params)?;

        let pj_ellps = if let Some(radius) = params.get("R") {
            Ellipsoid::sphere(radius.try_into()?)
        } else {
            // Retrieve datum
            let (method, ellps) = match params.get("datum") {
                Some(datum) => {
                    let datum =
                        datums::find_datum(datum.try_into()?).ok_or(Error::InvalidDatumError)?;
                    (Some(&datum.shift), Some(datum.ellps))
                }
                None => (None, None),
            };

            // Override Ellipsoid ?
            let ellps = match params.get("ellps") {
                Some(ellps) => ellipsoids::find_ellipsoid(ellps.try_into()?)
                    .ok_or(Error::InvalidEllipsoidError)?,
                None => ellps.unwrap_or(&ellipsoids::constants::WGS84),
            };

            Ellipsoid::from_ellipsoid_with_params(ellps, &params)
        }?;

        // Datum

        Ok(Self { pj_pm, pj_ellps })
    }

    /// Create projection from string
    pub fn from_str(s: &str) -> Result<Self> {
        Self::try_new(projstring::parse(s)?)
    }
}
