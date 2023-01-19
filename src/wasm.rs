//!
//! Wasm bindgen entry point
//!
use crate::{errors, proj, transform};
use wasm_bindgen::prelude::*;

/// ----------------------------
/// Wrapper for Projection
/// ---------------------------
#[wasm_bindgen]
pub struct Projection {
    inner: proj::Projection,
}

#[wasm_bindgen]
impl Projection {
    #[wasm_bindgen(constructor)]
    pub fn new(defn: &str) -> Result<Projection, JsError> {
        Ok(Self {
            inner: proj::Projection::from_projstr(defn)?,
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
    transform(src, dst, point)
}
