//! Implemented converions
//!
//! Reference: <https://proj.org/en/9.3/operations/conversions>
//!
//! Example:
//! ```
//! use proj4rs::conversions::Conversion;
//!
//! let conversion = Conversion::from_proj_string("+proj=axisswap +order=2,1").unwrap();
//! let mut points = (1., 2., 0.);
//! conversion.convert(&mut points).unwrap();
//! assert_eq!((2., 1., 0.), points);
//! ```

mod conversion;
pub use conversion::Conversion;

mod convert;
pub(crate) use convert::{Convert, ConvertParameters};

mod axisswap;
pub use axisswap::{AxisswapConversion, AxisswapOrdering};

mod noop;
pub use noop::NoopConversion;
