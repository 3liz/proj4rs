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
