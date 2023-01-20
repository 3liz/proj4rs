//!
//! Coordinate transformation library
//!
//! Based on Proj4js port of Proj4
//!
//! References:
//! * <http://docs.opengeospatial.org/as/18-005r5/18-005r5.html>
//! * <https://proj.org/development/reference/cpp/cpp_general.html>
//!
mod consts;
mod datum_params;
mod datum_transform;
mod datums;
mod ellipsoids;
mod ellps;
mod geocent;
mod nadgrids;
mod parameters;
mod prime_meridians;
mod projections;
mod projstring;
mod units;
mod utils;

pub mod adaptors;
pub mod errors;
pub mod proj;
pub mod transform;

// Include wasm entry point for wasm32-unknown-unknown
#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;

#[cfg(test)]
mod tests;
