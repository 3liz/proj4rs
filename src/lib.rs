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
//! The goal of proj4rs is not to be a remplacement of proj, but instead being a light
//! weight implementation of transformations from crs to crs that could be used
//! in WASM environment
//!
//! ## Usage
//!
//! Note that angular units are in radians, not degrees !
//!
//! Radian is natural unit for trigonometric operations, like proj, proj4rs use radians
//! for its operation while degrees are mostly used as end user input/output.
//!
//! Example:
//! ```
//! use proj4rs::proj::Proj;
//!
//! // EPSG:5174 - Example
//! let from = Proj::from_proj_string(concat!(
//!    "+proj=tmerc +lat_0=38 +lon_0=127.002890277778",
//!    " +k=1 +x_0=200000 +y_0=500000 +ellps=bessel",
//!    " +towgs84=-145.907,505.034,685.756,-1.162,2.347,1.592,6.342",
//!    " +units=m +no_defs +type=crs"
//! ))
//! .unwrap();
//!
//! // EPSG:4326 - WGS84, known to us as basic longitude and latitude.
//! let to = Proj::from_proj_string(concat!(
//!    "+proj=longlat +ellps=WGS84",
//!    " +datum=WGS84 +no_defs"
//! ))
//! .unwrap();
//!
//! let mut point_3d = (198236.3200000003, 453407.8560000006, 0.0);
//! proj4rs::transform::transform(&from, &to, &mut point_3d).unwrap();
//!
//! // XXX Note that angular unit is radians, not degrees !
//! point_3d.0 = point_3d.0.to_degrees();
//! point_3d.1 = point_3d.1.to_degrees();
//!
//! // Output in longitude, latitude, and height.
//! println!("{} {}",point_3d.0, point_3d.1); // 126.98069676435814, 37.58308534678718
//! ```
//!
//! ## Optional features
//!
//! * **geo-types**: [geo-types](<https://docs.rs/geo-types/latest/geo_types/>) support
//! * **logging**: support for logging with [log](https://docs.rs/log/latest/log/) crate.
//!   If activated for WASM, it will use the [console-log](https://docs.rs/console_log/latest/console_log/)
//!   adaptor.
//! * **wasm-strict**: used with WASM; Transformation operation will return exception as soon as we
//! have invalid coordinates or that the reprojection failed.
//! The default is to use a relaxed-mode that return NaN in case of projection failure: this is expected
//!   mostly from js app (at least with OpenLayer).
//! * **multi-thread**: Support for multi-thread with NAD Grid processing, this is activated by
//!   default and disabled when compiling for WASM.
//!
//! ## WKT Support
//!
//! There is no actual default support for WKT in proj4rs
//! If you are looking for WTK/Proje string conversion support in Rust,
//! then have a look at:
//!
//! - <https://github.com/3liz/proj4wkt-rs>
//! - <https://github.com/frewsxcv/crs-definitions>
//!
//! Note that the proj library provides a great implementation of the standard.
//!
//! ## Grid shift supports
//!
//! Nadgrid support is still experimental.
//! Currently, only Ntv2 multi grids are supported for native build and WASM.
//!

mod datum_params;
mod datum_transform;
mod datums;
mod ellipsoids;
mod ellps;
mod geocent;
mod math;
mod parameters;
pub(crate) use parameters::ParamList;
mod parse;
mod prime_meridians;
mod projstring;
mod units;

pub mod conversions;
pub(crate) use conversions::*;

pub mod adaptors;
pub mod errors;
pub(crate) use errors::Error as ProjError;
pub(crate) use errors::Result as ProjResult;

pub mod nadgrids;
pub mod proj;
pub mod projections;
pub mod transform;
pub(crate) use transform::Transform;

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
