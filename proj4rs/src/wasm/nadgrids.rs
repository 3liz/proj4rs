//!
//! WASM provider for nadgrids
//!
//! Use JS Dataview for passing nadgrids definition
//!
use js_sys::Uint8Array;
use wasm_bindgen::prelude::*;
use std::io::Cursor;

use crate::nadgrids::{catalog, formats};

#[wasm_bindgen]
pub fn add_nadgrid(key: &str, array: &Uint8Array) -> Result<(), JsError> {
    // NOTE: Is there a way to get a reference to the inner slice
    // of the array without copying the data ?
    let mut buffer = Cursor::new(array.to_vec());
    catalog::with_catalog(|cat| formats::read(cat, key, &mut buffer))?;
    Ok(())
}
