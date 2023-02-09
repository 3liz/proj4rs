//!
//! Coordinate transformation library
//!
//! Based on Proj4js port of Proj4
//!
//! References:
//! * <http://docs.opengeospatial.org/as/18-005r5/18-005r5.html>
//! * <https://proj.org/development/reference/cpp/cpp_general.html>
//!

mod datum_params;
mod datum_transform;
mod datums;
mod ellipsoids;
mod ellps;
mod geocent;
mod math;
mod nadgrids;
mod parameters;
mod parse;
mod prime_meridians;
mod projections;
mod projstring;
mod units;

pub mod adaptors;
pub mod errors;
pub mod proj;
pub mod transform;

// Reexport
pub use parameters::{ParamList, Parameter};
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

    pub(crate) use {
        __debug__ as debug, __error__ as error, __info__ as info, __trace__ as trace,
        __warn__ as warn,
    };
}
