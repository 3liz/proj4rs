//!
//! Define constants and derived constants
//!
//! Note: most of the definitions and constants are the same
//! as those defined in proj4js.
//!

// Note that TAU is 2*PI
// see https://doc.rust-lang.org/std/f64/consts/constant.TAU.html
pub use std::f64::consts::{FRAC_PI_2, PI, TAU};

// Was defined in proj4js for preventing divergence
// of Mollweied algorithm
pub const EPSLN: f64 = 1.0e-10;
pub const SEC_TO_RAD: f64 = 4.84813681109535993589914102357e-6;
