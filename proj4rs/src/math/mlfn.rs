//!
//! mlfn
//!  Meridional distance
//!
//!
use crate::errors::{Error, Result};

//  XXX Use clenshaw coefficients
//  with the third flattening ?
//  (cf Proj 9)

/// Alias for mlfn coefficients
pub(crate) type Enfn = (f64, f64, f64, f64, f64);

/// Meridional distance for ellipsoid and inverse
/// 8th degree - accurate to < 1e-5 meters when used in conjunction
/// with typical major axis values.
/// Inverse determines phi to EPS (1e-11) radians, about 1e-6 seconds.
pub(crate) fn enfn(es: f64) -> Enfn {
    const C00: f64 = 1.;
    const C02: f64 = 0.25;
    const C04: f64 = 0.046875;
    const C06: f64 = 0.01953125;
    const C08: f64 = 0.01068115234375;
    const C22: f64 = 0.75;
    const C44: f64 = 0.46875;
    const C46: f64 = 0.013_020_833_333_333_334;
    const C48: f64 = 0.007_120_768_229_166_667;
    const C66: f64 = 0.364_583_333_333_333_3;
    const C68: f64 = 0.005_696_614_583_333_334;
    const C88: f64 = 0.3076171875;

    let t = es * es;
    (
        C00 - es * (C02 + es * (C04 + es * (C06 + es * C08))),
        es * (C22 - es * (C04 + es * (C06 + es * C08))),
        t * (C44 - es * (C46 + es * C48)),
        t * es * (C66 - es * C68),
        t * t * es * C88,
    )
}

pub(crate) fn mlfn(phi: f64, mut sphi: f64, mut cphi: f64, en: Enfn) -> f64 {
    cphi *= sphi;
    sphi *= sphi;
    en.0 * phi - cphi * (en.1 + sphi * (en.2 + sphi * (en.3 + sphi * en.4)))
}

pub(crate) fn inv_mlfn(arg: f64, es: f64, en: Enfn) -> Result<f64> {
    const MAX_ITER: usize = 10;
    const EPS: f64 = 1e-11;
    let k = 1. / (1. - es);
    let mut phi = arg;
    let mut i = MAX_ITER;
    // rarely goes over 2 iterations
    while i > 0 {
        let s = phi.sin();
        let mut t = 1. - es * s * s;
        t = (mlfn(phi, s, phi.cos(), en) - arg) * (t * t.sqrt()) * k;
        phi -= t;
        if t.abs() < EPS {
            break;
        }
        i -= 1;
    }
    if i > 0 {
        Ok(phi)
    } else {
        Err(Error::InvMeridDistConvError)
    }
}
