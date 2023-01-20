//!
//! Define constants and derived constants
//!
//! Note: most of the definitions and constants are the same
//! as those defined in proj4js.
//!

// Note that TAU is 2*PI
// see https://doc.rust-lang.org/std/f64/consts/constant.TAU.html
pub(crate) use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

// Was defined in proj4js for preventing divergence
// of Mollweied algorithm
pub(crate) const EPS_10: f64 = 1.0e-10;

// Other value op epsilon used
pub(crate) const EPS_12: f64 = 1.0e-12;

// XXX float has excessive precision
//pub const SEC_TO_RAD: f64 = 4.84813681109535993589914102357e-6;
pub(crate) const SEC_TO_RAD: f64 = 4.848_136_811_095_36e-6;
