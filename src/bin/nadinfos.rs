//!
//! Display infos about NAD grid
//!
use proj4rs::errors::{Error, Result};
use proj4rs::nadgrids::{files::read_from_file, Catalog};
use std::env;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 {
        println!("Usage: nadinfos <path>");
        return Err(Error::InvalidParameterValue("Missing filename"));
    }

    let key = &args[1];

    let catalog = Catalog::default();
    read_from_file(&catalog, key)?;

    match catalog.find(key) {
        Some(iter) => {
            iter.for_each(|g| println!("{g}#"));
            Ok(())
        }
        None => Err(Error::NadGridNotAvailable),
    }
}
