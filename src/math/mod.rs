//!
//! Utilities
//!
//!
use crate::errors::{Error, Result};

pub(crate) mod consts {
    //!
    //! Define constants
    //!

    // Note that TAU is 2*PI
    // see https://doc.rust-lang.org/std/f64/consts/constant.TAU.html
    pub(crate) use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

    // Was defined in proj4js for preventing divergence
    // of Mollweied algorithm
    pub(crate) const EPS_10: f64 = 1.0e-10;

    // Other value op epsilon used
    pub(crate) const EPS_12: f64 = 1.0e-12;

    // Other value op epsilon used
    pub(crate) const EPS_7: f64 = 1.0e-7;

    // XXX float has excessive precision
    //pub const SEC_TO_RAD: f64 = 4.84813681109535993589914102357e-6;
    pub(crate) const SEC_TO_RAD: f64 = 4.848_136_811_095_36e-6;
}

// Redefinition of mathematical functions
//
// Some of these functions has been redefined for various reason.
// It would be nice to investigate if some of them are still relevant
//
// Note that proj redefine ln1p (i.e ln(1+x)), while rust rely on platform native (libm)
// implementation:
//
// ```C
// static double log1py(double x) {              /* Compute log(1+x) accurately */
//    volatile double
//      y = 1 + x,
//      z = y - 1;
//    /* Here's the explanation for this magic: y = 1 + z, exactly, and z
//     * approx x, thus log(y)/z (which is nearly constant near z = 0) returns
//     * a good approximation to the true log(1 + x)/x.  The multiplication x *
//    * (log(y)/z) introduces little additional error. */
//    return z == 0 ? x : x * log(y) / z;
// ```
//
// For now we are going to stick to the native implementation of `ln_1p`, let's see if that
// may cause problems in the future
//
//
// The same for hypot, for now we are going to stick to the native implementation.
// since latest version of glibc seems to handle case of potential overflow.

//  ----------
//  asinh
//  ---------
//
// In the case of [`asinh`], rust define this as (https://doc.rust-lang.org/src/std/f64.rs.html#882-884)
//
// ```rust
// pub fn asinh(self) -> f64 {
//        (self.abs() + ((self * self) + 1.0).sqrt()).ln().copysign(self)
//  }
// ```
//
// Note that proj use the following formula:
// ```rust
// #[inline]
// pub fn asinh(x: f64) -> f64 {
//     let y = x.abs();         // Enforce odd parity
//     (y * (1. + y/(hypot(1.0,y) + 1.))).ln_1p().copysign(x)
// }
// ```

// The formula below is mathematically equivalent, but the rust version use
// a naive implementation of `hypot`
// wich may (eventually) leads to overflow.
//
// We prefer to use our own implementation using [`hypot`] with the simpler
// rust formula. This implementation will give accurate result for `0.89e308f64` while the
// `[f64::asinh`] implementation overflow (return `f64::INFINITE`)
#[inline]
pub fn asinh(x: f64) -> f64 {
    (x.abs() + 1.0f64.hypot(x)).ln().copysign(x)
}

//pub fn asinh(x: f64) -> f64 {
//    let y = x.abs();         // Enforce odd parity
//    (y * (1. + y/(1.0f64.hypot(y) + 1.))).ln_1p().copysign(x)
//}

mod adjlon;
mod gauss;
mod mlfn;
mod msfn;
mod phi2;
mod qsfn;
mod tsfn;

pub(crate) use adjlon::adjlon;
pub(crate) use gauss::{gauss, gauss_ini, inv_gauss, Gauss};
pub(crate) use mlfn::{enfn, inv_mlfn, mlfn, Enfn};
pub(crate) use msfn::msfn;
pub(crate) use phi2::phi2;
pub(crate) use qsfn::qsfn;
pub(crate) use tsfn::tsfn;
