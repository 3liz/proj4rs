//!
//! Coordinate transformation library
//!
//! Based on Proj4js port of Proj4
//!
//! References:
//! * <http://docs.opengeospatial.org/as/18-005r5/18-005r5.html>
//! * <https://proj.org/development/reference/cpp/cpp_general.html>
//!
pub(crate) mod constants;
pub(crate) mod datum_params;
pub(crate) mod datum_transform;
pub(crate) mod datums;
pub(crate) mod ellipsoids;
pub(crate) mod ellps;
pub(crate) mod geocent;
pub(crate) mod nadgrids;
pub(crate) mod parameters;
pub(crate) mod prime_meridians;
pub(crate) mod projstring;
pub(crate) mod transform;
pub(crate) mod units;

pub mod errors;
pub mod proj;

// Include wasm entry point for wasm32-unknown-unknown
#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm;
