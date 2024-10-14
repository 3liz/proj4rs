//!
//! Wasm bindgen entry point
//!
mod nadgrids;

use crate::{errors, proj, transform, transform::TransformClosure};
use wasm_bindgen::prelude::*;

use crate::log;

// Js entry point
#[cfg(feature = "with-wasm-entrypoint")]
#[wasm_bindgen(start)]
pub fn main() {
    #[cfg(feature = "logging")]
    console_log::init_with_level(log::Level::Trace).unwrap();

    log::info!("Initialized proj4rs wasm module.")
}

// ----------------------------
// Wrapper for Projection
// ---------------------------
#[wasm_bindgen]
pub struct Projection {
    inner: proj::Proj,
}

#[wasm_bindgen]
impl Projection {
    #[wasm_bindgen(constructor)]
    pub fn new(defn: &str) -> Result<Projection, JsError> {
        Ok(Self {
            inner: proj::Proj::from_user_string(defn)?,
        })
    }

    #[wasm_bindgen(getter, js_name = projName)]
    pub fn projname(&self) -> String {
        self.inner.projname().into()
    }

    #[wasm_bindgen(getter, js_name = isLatlon)]
    pub fn is_latlong(&self) -> bool {
        self.inner.is_latlong()
    }

    #[wasm_bindgen(getter, js_name = isGeocentric)]
    pub fn is_geocent(&self) -> bool {
        self.inner.is_geocent()
    }

    #[wasm_bindgen(getter)]
    pub fn axis(&self) -> String {
        String::from_utf8_lossy(self.inner.axis()).into_owned()
    }

    #[wasm_bindgen(getter, js_name = isNormalizedAxis)]
    pub fn is_normalized_axis(&self) -> bool {
        self.inner.is_normalized_axis()
    }

    #[wasm_bindgen(getter)]
    pub fn to_meter(&self) -> f64 {
        self.inner.to_meter()
    }

    #[wasm_bindgen(getter)]
    pub fn units(&self) -> String {
        self.inner.units().into()
    }
}

// ----------------------------
// Wrapper for Transform
// ---------------------------
#[wasm_bindgen]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[wasm_bindgen]
impl Point {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Self { x, y, z }
    }
}

impl transform::Transform for Point {
    /// Strict mode: return exception
    /// as soon as with have invalid coordinates or
    /// that the reprojection failed
    #[cfg(feature = "wasm-strict")]
    fn transform_coordinates<F: TransformClosure>(&mut self, f: &mut F) -> errors::Result<()> {
        f(self.x, self.y, self.z).map(|(x, y, z)| {
            self.x = x;
            self.y = y;
            self.z = z;
        })
    }
    /// Relaxed mode: allow transformation failure: return NAN in case
    /// of projection failure
    /// Note: this is what is expected mostly from js app (at least OpenLayer)
    #[cfg(not(feature = "wasm-strict"))]
    fn transform_coordinates<F: TransformClosure>(&mut self, f: &mut F) -> errors::Result<()> {
        f(self.x, self.y, self.z)
            .map(|(x, y, z)| {
                self.x = x;
                self.y = y;
                self.z = z;
            })
            .or_else(|_err| {
                // This will be activated with 'logging' feature
                log::error!("{:?}: ({}, {}, {})", _err, self.x, self.y, self.z);
                self.x = f64::NAN;
                self.y = f64::NAN;
                self.z = f64::NAN;
                Ok(())
            })
    }
}

#[wasm_bindgen]
pub fn transform(src: &Projection, dst: &Projection, point: &mut Point) -> Result<(), JsError> {
    // log::debug!("transform: {}, {}, {}", point.x, point.y, point.z);

    if point.x.is_nan() || point.y.is_nan() {
        return Err(JsError::from(errors::Error::NanCoordinateValue));
    }

    if src.inner.is_latlong() {
        point.x = point.x.to_radians();
        point.y = point.y.to_radians();
    }
    transform::transform(&src.inner, &dst.inner, point)?;
    if dst.inner.is_latlong() {
        point.x = point.x.to_degrees();
        point.y = point.y.to_degrees();
    }
    Ok(())
}
