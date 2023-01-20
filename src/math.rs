//!
//! Utilities
//!
//!
use crate::consts::{EPS_12, FRAC_PI_2, PI, TAU};

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
    (0.5 * (FRAC_PI_2 - phi)).tan() / ((1. - sinphi * e) / (1. + sinphi * e)).powf(0.5 * e)
}
