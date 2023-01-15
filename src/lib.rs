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
pub(crate) mod datums;
pub(crate) mod ellipsoids;
pub(crate) mod ellps;
pub(crate) mod parameters;
pub(crate) mod prime_meridians;
pub(crate) mod projstring;
pub(crate) mod units;

pub mod errors;
pub mod projection;
