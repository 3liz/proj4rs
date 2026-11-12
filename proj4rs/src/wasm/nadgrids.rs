//!
//! WASM provider for nadgrids
//!
//! Use JS Dataview for passing nadgrids definition
//!
use js_sys::Uint8Array;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use wasm_bindgen::prelude::*;

use crate::nadgrids::{catalog, formats};

/// Buffer size used for reading from the JS array
const READ_BUFFER_SIZE: usize = 64 * 1024;

/// Streaming reader over a JS `Uint8Array`
///
/// JS memory cannot be borrowed from the wasm linear memory, so data has to
/// be copied at least once. This reader copies chunks directly from the JS
/// buffer into the destination buffer, avoiding an intermediate copy of
/// the whole array in wasm memory.
struct JsArrayReader<'a> {
    array: &'a Uint8Array,
    pos: u32,
    len: u32,
}

impl<'a> JsArrayReader<'a> {
    fn new(array: &'a Uint8Array) -> Self {
        Self {
            array,
            pos: 0,
            len: array.length(),
        }
    }
}

impl Read for JsArrayReader<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = (self.len - self.pos).min(u32::try_from(buf.len()).unwrap_or(u32::MAX));
        if n > 0 {
            self.array
                .subarray(self.pos, self.pos + n)
                .copy_to(&mut buf[..n as usize]);
            self.pos += n;
        }
        Ok(n as usize)
    }
}

impl Seek for JsArrayReader<'_> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let new = match pos {
            SeekFrom::Start(p) => i64::try_from(p).unwrap_or(i64::MAX),
            SeekFrom::End(o) => (self.len as i64).saturating_add(o),
            SeekFrom::Current(o) => (self.pos as i64).saturating_add(o),
        };
        if new < 0 || new > self.len as i64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid seek position",
            ));
        }
        self.pos = new as u32;
        Ok(self.pos as u64)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.pos as u64)
    }
}

#[wasm_bindgen]
pub fn add_nadgrid(key: &str, array: &Uint8Array) -> Result<(), JsError> {
    // Each read crosses the JS boundary: buffer it so that the many small
    // reads from the parsers are turned into a few large ones.
    let mut reader = BufReader::with_capacity(READ_BUFFER_SIZE, JsArrayReader::new(array));
    catalog::with_catalog(|cat| formats::read(cat, key, &mut reader))?;
    Ok(())
}
