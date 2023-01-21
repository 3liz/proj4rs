//!
//! Utilities
//!
//!
use crate::consts::{EPS_10, EPS_12, FRAC_PI_2, PI, TAU};
use crate::errors::{Error, Result};

pub(crate) fn adjlon(mut lon: f64) -> f64 {
    // Let lon slightly overshoot,
    // to avoid spurious sign switching at the date line
    if lon.abs() >= PI + EPS_12 {
        // adjust to 0..2pi rad
        lon += PI;

        // remove integral # of 'revolutions'
        lon -= TAU * (lon / TAU).floor();

        // adjust back to -pi..pi rad
        lon -= PI;
    }
    lon
}

#[inline(always)]
pub(crate) fn msfn(sinphi: f64, cosphi: f64, es: f64) -> f64 {
    cosphi / (1. - es * sinphi * sinphi).sqrt()
}

#[inline(always)]
pub(crate) fn tsfn(phi: f64, sinphi: f64, e: f64) -> f64 {
    //  XXX Avoid division by zero, check denominator
    (0.5 * (FRAC_PI_2 - phi)).tan() / ((1. - sinphi * e) / (1. + sinphi * e)).powf(0.5 * e)
}

const PHI2_NITER: i32 = 15;
/// Determine latitude angle phi-2.
/// Inputs:
///  ts = exp(-psi) where psi is the isometric latitude (dimensionless)
///  e = eccentricity of the ellipsoid (dimensionless)
/// Output:
/// phi = geographic latitude (radians)
/// Here isometric latitude is defined by
/// psi = log( tan(pi/4 + phi/2) *
///            ( (1 - e*sin(phi)) / (1 + e*sin(phi)) )^(e/2) )
///      = asinh(tan(phi)) - e * atanh(e * sin(phi))
/// This routine inverts this relation using the iterative scheme given
/// by Snyder (1987), Eqs. (7-9) - (7-11)
///
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
