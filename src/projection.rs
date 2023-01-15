//!
//! Derived values from projection definition.
//!
//! Precomputed for optimization
//!

use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::{datums, ellipsoids, projstring};

pub struct Projection {
    pj_ellps: Ellipsoid,
}

impl Projection {
    /// Consume a  ParamList and create a Projection object
    pub fn try_new(params: ParamList) -> Result<Self> {
        // Projection name
        let projname = params.get("proj").ok_or(Error::MissingProjectionError);

        // R takes precedence (from proj, not in proj4js)
        let pj_ellps = if let Some(radius) = params.get("R") {
            Ellipsoid::sphere(radius.try_into()?)
        } else {
            // Retrieve datum
            let (method, ellps) = match params.get("datum") {
                Some(datum) => {
                    let (method, ellps) =
                        datums::datum_defn(datum.try_into()?).ok_or(Error::InvalidDatumError)?;
                    (Some(method), Some(ellps))
                }
                None => (None, None),
            };

            // Override Ellipsoid ?
            let ellps = match params.get("ellps") {
                Some(ellps) => {
                    ellipsoids::ellps_defn(ellps.try_into()?).ok_or(Error::InvalidEllipsoidError)?
                }
                None => ellps.unwrap_or(&ellipsoids::constants::WGS84),
            };

            Ellipsoid::from_ellipsoid_with_params(ellps, &params)
        }?;

        Ok(Self { pj_ellps })
    }

    /// Create projection from string
    pub fn from_str(s: &str) -> Result<Self> {
        Self::try_new(projstring::parse(s)?)
    }
}
