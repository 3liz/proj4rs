//!
//! Read grid from files
//!
use std::env;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::errors::{Error, Result};
use crate::nadgrids::Catalog;

mod ntv2;

use super::header::Header;
use ntv2::read_ntv2;

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

pub(crate) enum FileType {
    Ntv1,
    Ntv2,
    Gtx,
    Ctable,
    Ctable2,
}

/// Recognize grid file type
pub(crate) fn recognize<R: Read + Seek>(key: &str, read: &mut R) -> Result<FileType> {
    const BUFSIZE: usize = 160;
    let pos = read.stream_position()?;
    let mut header = Header::<BUFSIZE>::new();

    let rv = header.read_partial(read).map(|size| {
        if size >= 144 + 16
            && header.cmp_str(0, "HEADER")
            && header.cmp_str(96, "W_GRID")
            && header.cmp_str(144, "TO      NAD83   ")
        {
            FileType::Ntv1
        } else if size >= 48 + 7 && header.cmp_str(0, "NUM_OREC") && header.cmp_str(48, "GS_TYPE") {
            FileType::Ntv2
        } else if key.ends_with("gtx") || key.ends_with("GTX") {
            FileType::Gtx
        } else if size >= 9 && header.cmp_str(0, "CTABLE V2") {
            FileType::Ctable2
        } else {
            // Ctable fallback
            FileType::Ctable
        }
    });

    // Restore position
    read.seek(SeekFrom::Start(pos))?;
    rv
}

/// Grid builder that read from a file
pub fn read_from_file(catalog: &Catalog, key: &str) -> Result<()> {
    // Use a BufReader for efficiency
    read(
        catalog,
        key,
        &mut BufReader::new(File::open(default_file_finder(key)?)?),
    )
}

/// Read a grid from a file given by `key`
pub fn read<R: Read + Seek>(catalog: &Catalog, key: &str, read: &mut R) -> Result<()> {
    // Guess the file
    match recognize(key, read)? {
        FileType::Ntv2 => read_ntv2(catalog, key, read),
        _ => Err(Error::UnknownGridFormat),
    }
}
