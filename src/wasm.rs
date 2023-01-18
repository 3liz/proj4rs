//!
//! Wasm bindgen entry point
//!
use crate::errors::Error;
use crate::proj;
use wasm_bindgen::prelude::*;

/// Wrapper for Projection
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
