//!
//! Gauss
//!

// Original copyright
//
// Copyright (c) 2003   Gerald I. Evenden
//
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
use super::consts::{FRAC_PI_2, FRAC_PI_4};
use crate::errors::{Error, Result};

#[inline]
fn srat(esinp: f64, ratexp: f64) -> f64 {
    ((1. - esinp) / (1. + esinp)).powf(ratexp)
}

#[derive(Debug)]
pub(crate) struct Gauss {
    c: f64,
    k: f64,
    e: f64,
    ratexp: f64,
}

pub(crate) fn gauss_ini(e: f64, phi0: f64) -> Result<(Gauss, f64, f64)> {
    let es = e * e;
    let (sphi, mut cphi) = phi0.sin_cos();

    cphi *= cphi;

    let rc = (1. - es).sqrt() / (1. - es * sphi * sphi);
    let c = (1. + es * cphi * cphi / (1. - es)).sqrt();

    if c == 0. {
        return Err(Error::ToleranceConditionError);
    }

    let chi = (sphi / c).asin();
    let ratexp = 0.5 * c * e;
    let k = (0.5 * chi + FRAC_PI_4).tan()
        / ((0.5 * phi0 + FRAC_PI_4).tan().powf(c) * srat(e * sphi, ratexp));
    Ok((Gauss { c, k, e, ratexp }, chi, rc))
}

pub(crate) fn gauss(lam: f64, phi: f64, en: &Gauss) -> (f64, f64) {
    (
        // lam
        en.c * lam,
        // phi
        2. * (en.k * (0.5 * phi + FRAC_PI_4).tan().powf(en.c) * srat(en.e * phi.sin(), en.ratexp))
            .atan()
            - FRAC_PI_2,
    )
}

pub(crate) fn inv_gauss(lam: f64, mut phi: f64, en: &Gauss) -> Result<(f64, f64)> {
    const DEL_TOL: f64 = 1.0e-14;
    const MAX_ITER: usize = 20;
    let mut i = MAX_ITER;
    let num = ((0.5 * phi + FRAC_PI_4).tan() / en.k).powf(1. / en.c);
    // XXX should try
    /*
    match (0..MAX_ITER).try_fold(phi, |phi, _| {
        let e_phi = 2. * (num * srat(en.e * phi.sin(), -0.5 * en.e)).atan() - FRAC_PI_2;
        if (e_phi - phi).abs() < DEL_TOL {
            Break(e_phi)
        }
        Continue(e_phi)
    }) {
        Break(phi) => Ok((lam / en.c, phi)),
        _ => Err(Error::InvMeridDistConvError),
    }
    */
    while i > 0 {
        let e_phi = 2. * (num * srat(en.e * phi.sin(), -0.5 * en.e)).atan() - FRAC_PI_2;
        if (e_phi - phi).abs() < DEL_TOL {
            phi = e_phi;
            break;
        }
        phi = e_phi;
        i -= 1;
    }
    if i > 0 {
        Ok((lam / en.c, phi))
    } else {
        Err(Error::InvMeridDistConvError)
    }
}
