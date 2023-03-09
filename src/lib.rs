//!
//! Coordinate transformation library
//!
//! Based on Proj4 implementation
//!
//! References:
//! * <http://docs.opengeospatial.org/as/18-005r5/18-005r5.html>
//! * <https://proj.org/development/reference/cpp/cpp_general.html>
//!
//! The aim of Proj4rs is to provide at short term the same functionality as the original
//! proj4 library.
//!
//! The long term project is to integrate feature from the proj library in its latest
//! version.
//!
//! The goal of proj4rs is not to be a remplacement of proj, but instead beeing a light
//! weight implementation of transformations from crs to crs that could be used
//! in WASM environment
//!
//! There is no actual support for WKT, if such supports would exist one day, it would be under
//! a dedicated crate for transforming proj string to to WKT and vice-versa.
//!
//! If you need full support for WKT, please rely on proj which provide
//! a great implementation of the standard.
//!

mod datum_params;
mod datum_transform;
mod datums;
mod ellipsoids;
mod ellps;
mod geocent;
mod math;
mod parameters;
mod parse;
mod prime_meridians;
mod projstring;
mod units;

pub mod adaptors;
pub mod errors;
pub mod nadgrids;
pub mod proj;
pub mod projections;
pub mod transform;

// Reexport
pub use proj::Proj;

// Include wasm entry point for wasm32-unknown-unknown
#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;

#[cfg(test)]
mod tests;

#[cfg(target_arch = "wasm32")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

// log for logging (optional).
#[cfg(feature = "logging")]
use log;

#[cfg(not(feature = "logging"))]
mod log {
    // Use __XXX__ to prevent 'ambiguous name' error
    // when exporting
    macro_rules! __trace__    ( ($($tt:tt)*) => {{}} );
    macro_rules! __debug__    ( ($($tt:tt)*) => {{}} );
    macro_rules! __error__    ( ($($tt:tt)*) => {{}} );
    macro_rules! __info__     ( ($($tt:tt)*) => {{}} );
    macro_rules! __warn__     ( ($($tt:tt)*) => {{}} );

    #[allow(unused_imports)]
    pub(crate) use {
        __debug__ as debug, __error__ as error, __info__ as info, __trace__ as trace,
        __warn__ as warn,
    };
}
