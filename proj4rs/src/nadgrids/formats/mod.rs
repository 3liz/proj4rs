//!
//! Read grid from files
//!
use std::io::{Read, Seek, SeekFrom};

use crate::errors::{Error, Result};
use crate::nadgrids::Catalog;

#[cfg(feature = "tiff")]
mod geotiff;
mod ntv2;

#[cfg(feature = "tiff")]
use geotiff::is_tiff;

// Dummy tiff check
#[cfg(not(feature = "tiff"))]
fn is_tiff<const N: usize>(_: usize, _: Header<N>) -> bool {
    false
}

use super::header::Header;

pub(crate) enum FileType {
    Ntv1,
    Ntv2,
    Gtx,
    Ctable2,
    Tiff,
    Unknown,
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
        } else if is_tiff(size, header) {
            FileType::Tiff
        } else {
            FileType::Unknown
        }
    });

    // Restore position
    read.seek(SeekFrom::Start(pos))?;
    rv
}

/// Read a grid from IO stream
pub fn read<R: Read + Seek>(catalog: &Catalog, key: &str, read: &mut R) -> Result<()> {
    // Guess the file
    match recognize(key, read)? {
        FileType::Ntv2 => ntv2::read_ntv2(catalog, key, read),
        #[cfg(feature = "tiff")]
        FileType::Tiff => geotiff::read_tiff(catalog, key, read),
        _ => Err(Error::UnknownGridFormat),
    }
}
