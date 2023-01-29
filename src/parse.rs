//!
//! Parser for numbers
//!
//! If we are in wasm mode, then fallback to the Js functions
//! this is a 20Ko gain for the .wasm binary.
//!

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
mod wasm {
    //! Use js lib

    use crate::errors::Error;
    use js_sys::{parse_float, parse_int};

    pub trait FromStr: Sized {
        type Err;

        fn from_str(s: &str) -> Result<Self, Self::Err>;
    }

    impl FromStr for f64 {
        type Err = Error;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let v = parse_float(s);
            if v.is_nan() {
                Err(Error::JsParseError)
            } else {
                Ok(v)
            }
        }
    }

    impl FromStr for i32 {
        type Err = Error;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let v = parse_int(s, 10);
            if v.is_nan() {
                Err(Error::JsParseError)
            } else {
                Ok(v as i32)
            }
        }
    }

    impl FromStr for bool {
        type Err = Error;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s {
                "true" => Ok(true),
                "false" => Ok(false),
                _ => Err(Error::JsParseError),
            }
        }
    }
}

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
pub use wasm::FromStr;

#[cfg(all(not(target_arch = "wasm32"), not(target_os = "unknown")))]
pub use std::str::FromStr;
