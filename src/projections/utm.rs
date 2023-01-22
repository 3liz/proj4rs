//!
//! Universal transverse Mercator UTM
//!
//! based on etmerc
//!
//! Parameters:
//!
//! proj: utm
//!
use crate::errors::{Error, Result};
use crate::math::{adjlon, consts::PI};
use crate::parameters::ParamList;
use crate::proj::ProjData;

use super::etmerc;

// Projection stub
super::projection!(utm, "utm");

#[derive(Debug)]
pub(crate) struct Projection(etmerc::Projection);

impl Projection {
    pub fn init(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        if p.lam0 < -1000. || p.lam0 > 1000. {
            return Err(Error::InvalidUtmZone);
        }

        p.x0 = 500_000.;
        p.y0 = if params.check_option("south")? {
            10_000_000.
        } else {
            0.
        };

        let zone = params.try_value::<u8>("zone").and_then(|zone| match zone {
            Some(zone) => {
                if (1..=60).contains(&zone) {
                    Ok(zone as f64)
                } else {
                    Err(Error::InvalidUtmZone)
                }
            }
            None => {
                // nearest central meridian input
                let zone = ((adjlon(p.lam0) + PI) * 30. / PI).floor().round();
                if (1. ..=60.).contains(&zone) {
                    Ok(zone)
                } else {
                    Err(Error::InvalidUtmZone)
                }
            }
        })?;

        p.lam0 = (zone + 0.5) * PI / 30. - PI;
        p.k0 = 0.9996;
        p.phi0 = 0.;

        Ok(Self(etmerc::Projection::init(p, params)?))
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        self.0.forward(lam, phi, z)
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        self.0.inverse(x, y, z)
    }

    pub const fn has_inverse() -> bool {
        etmerc::Projection::has_inverse()
    }

    pub const fn has_forward() -> bool {
        etmerc::Projection::has_forward()
    }
}
