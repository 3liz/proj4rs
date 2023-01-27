//! Determine latitude angle phi-2.
//!
//! Inputs:
//!  ts = exp(-psi) where psi is the isometric latitude (dimensionless)
//!  e = eccentricity of the ellipsoid (dimensionless)
//! Output:
//! phi = geographic latitude (radians)
//! Here isometric latitude is defined by
//! psi = log( tan(pi/4 + phi/2) *
//!            ( (1 - e*sin(phi)) / (1 + e*sin(phi)) )^(e/2) )
//!      = asinh(tan(phi)) - e * atanh(e * sin(phi))
//! This routine inverts this relation using the iterative scheme given
//! by Snyder (1987), Eqs. (7-9) - (7-11)
//!
use super::consts::{EPS_10, FRAC_PI_2};
use crate::errors::{Error, Result};

const PHI2_NITER: i32 = 15;

pub(crate) fn phi2(ts: f64, e: f64) -> Result<f64> {
    let eccnth = 0.5 * e;
    let mut phi = FRAC_PI_2 - 2. * ts.atan();
    let mut i = PHI2_NITER;
    while i > 0 {
        let con = e * phi.sin();
        let dphi = FRAC_PI_2 - 2. * (ts * ((1. - con) / (1. + con)).powf(eccnth)).atan() - phi;

        phi += dphi;

        if dphi.abs() <= EPS_10 {
            break;
        }

        i -= 1;
    }

    if i <= 0 {
        Err(Error::NonInvPhi2Convergence)
    } else {
        Ok(phi)
    }
}
