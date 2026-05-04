//!
//! Read grid from files
//!
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::errors::{Error, Result};
use crate::nadgrids::{Catalog, formats};

/// Define a default file finder functions
fn default_file_finder(name: &str) -> Result<PathBuf> {
    let p = Path::new(name);
    match p.exists().then_some(p.into()).or_else(|| {
        if let Ok(val) = env::var("PROJ_NADGRIDS").or_else(|_| env::var("PROJ_DATA")) {
            val.split(':').find_map(|s| {
                let p = Path::new(s).join(name);
                p.exists().then_some(p)
            })
        } else {
            None
        }
    }) {
        Some(p) => Ok(p),
        None => Err(Error::GridFileNotFound(name.into())),
    }
}

/// Grid builder that read from a file
pub fn read_from_file(catalog: &Catalog, key: &str) -> Result<()> {
    // Use a BufReader for efficiency
    formats::read(
        catalog,
        key,
        &mut BufReader::new(File::open(default_file_finder(key)?)?),
    )
}
