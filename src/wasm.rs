//!
//! Wasm bindgen entry point
//!
use crate::{errors, proj, transform};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

/// ----------------------------
/// Wrapper for Projection
/// ---------------------------
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
}

/// ----------------------------
/// Wrapper for Transform
/// ---------------------------
#[wasm_bindgen]
pub struct Point {
    x: f64,
    y: f64,
    z: f64,
}

#[wasm_bindgen]
impl Point {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Self { x, y, z }
    }

    #[wasm_bindgen(getter)]
    pub fn x(&self) -> f64 {
        self.x
    }
    #[wasm_bindgen(setter)]
    pub fn set_x(&mut self, x: f64) {
        self.x = x;
    }
    #[wasm_bindgen(getter)]
    pub fn y(&self) -> f64 {
        self.y
    }
    #[wasm_bindgen(setter)]
    pub fn set_y(&mut self, y: f64) {
        self.y = y;
    }
    #[wasm_bindgen(getter)]
    pub fn z(&self) -> f64 {
        self.z
    }
    #[wasm_bindgen(setter)]
    pub fn set_z(&mut self, z: f64) {
        self.z = z;
    }
}

impl transform::Transform for Point {
    fn transform_coordinates<F>(&mut self, mut f: F) -> errors::Result<()>
    where
        F: FnMut(f64, f64, f64) -> errors::Result<(f64, f64, f64)>,
    {
        f(self.x, self.y, self.z).map(|(x, y, z)| {
            self.x = x;
            self.y = y;
            self.z = z;
        })
    }
}

#[wasm_bindgen]
pub fn transform(src: &Projection, dst: &Projection, point: &mut Point) -> Result<(), JsError> {
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
