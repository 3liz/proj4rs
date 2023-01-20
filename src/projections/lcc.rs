//!
//! Lambert Conformal Conic
//!
//! Paramètres:
//!
//! proj: lcc
//!
//! lat_1:
//! lat_2:
//! lat_0:
//!

use crate::consts::{EPS_10, FRAC_PI_2, FRAC_PI_4};
use crate::errors::{Error, Result};
use crate::math::{msfn, tsfn};
use crate::parameters::ParamList;
use crate::proj::Proj;

use super::{ProjParams, ProjSetup};

pub(super) const NAME: &str = "lcc";

#[derive(Debug, Default)]
pub(crate) struct Projection {
    phi1: f64,
    phi2: f64,
    n: f64,
    rho0: f64,
    c: f64,
    ellips: bool,
}

impl Projection {
    pub fn init(p: &mut Proj, params: &ParamList) -> Result<ProjSetup> {
        let phi1 = params.try_value("lat_1")?.unwrap_or(0.);
        let phi2 = params.try_value("lat_2")?.unwrap_or_else(|| {
            p.phi0 = p.phi0.or(Some(phi1));
            phi1
        });

        if (phi1 + phi2).abs() < EPS_10 {
            return Err(Error::ProjErrConicLatEqual);
        }

        let phi0 = p.phi0();

        let sinphi = phi1.sin();
        let cosphi = phi1.cos();
        let secant = (phi1 - phi2).abs() >= EPS_10;

        let mut n = sinphi;
        let el = &p.ellps;

        let ellips = el.es != 0.;

        let (c, rho0);

        if ellips {
            let m1 = msfn(sinphi, cosphi, el.es);
            let ml1 = tsfn(phi1, sinphi, el.e);
            // secant zone
            if secant {
                let sinphi = phi2.sin();
                n = (m1 / msfn(sinphi, phi2.cos(), el.es)).ln();
                n /= (ml1 / tsfn(phi2, sinphi, el.e)).ln();
            }
            c = m1 * ml1.powf(-n) / n;
            rho0 = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
                0.
            } else {
                c * tsfn(phi0, phi0.sin(), el.e).powf(n)
            }
        } else {
            if secant {
                n = (cosphi / phi2.cos()).ln()
                    / ((FRAC_PI_4 + 0.5 * phi2).tan() / (FRAC_PI_4 + 0.5 * phi1).tan()).ln();
            }
            c = cosphi * (FRAC_PI_4 + 0.5 * phi1).tan().powf(n) / n;
            rho0 = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
                0.
            } else {
                c * (FRAC_PI_4 + 0.5 * phi0).tan().powf(-n)
            }
        }

        Ok((
            ProjParams::lcc(Self {
                phi1,
                phi2,
                n,
                rho0,
                c,
                ellips,
            }),
            Some(Self::inverse),
            Some(Self::forward),
        ))
    }

    pub fn forward(p: &Proj, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam / p.ellps.a, phi / p.ellps.a, z))
    }

    pub fn inverse(p: &Proj, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam * p.ellps.a, phi * p.ellps.a, z))
    }
}
