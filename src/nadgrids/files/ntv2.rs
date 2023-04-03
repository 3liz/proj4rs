//!
//! Nadgrid parser
//!
use crate::errors::{Error, Result};
use crate::log::trace;
use crate::math::consts::SEC_TO_RAD;
use crate::nadgrids::header::error_str::*;
use crate::nadgrids::header::{Endianness, Header};
use crate::nadgrids::{Catalog, Grid, GridId, Lp};
use std::io::Read;

const NTV2_HEADER_SIZE: usize = 11 * 16;

/// Ntv2 reader
pub(super) fn read_ntv2<R: Read>(catalog: &Catalog, key: &str, read: &mut R) -> Result<()> {
    let mut head = Header::<NTV2_HEADER_SIZE>::new();

    trace!("Reading  ntv2 {}", key);

    // Read overview header
    head.read(read)?;
    // Check endianness
    head.endian = if head.get_u8(8) == 11 {
        Endianness::native()
    } else {
        Endianness::other()
    };

    let nsubgrids = head.get_u32(40) as usize;

    trace!("Reading ntv2 {} subgrids {}", key, nsubgrids);

    // Read subsequent grids
    (0..nsubgrids).try_for_each(|_| read_ntv2_grid(catalog, key, head.read(read)?, read))
}

/// Read ntv2 grid data
fn read_ntv2_grid<R: Read>(
    catalog: &Catalog,
    key: &str,
    head: &Header<NTV2_HEADER_SIZE>,
    read: &mut R,
) -> Result<()> {
    match head.get_str(0, 8) {
        Ok("SUB_NAME") => Ok(()),
        _ => Err(Error::InvalidNtv2GridFormat(ERR_INVALID_HEADER)),
    }?;

    let id = head.get_id(8);
    let mut lineage = head.get_id(24);
    if lineage.as_str().trim() == "NONE" {
        lineage = GridId::root();
    }

    let mut ll = Lp {
        lam: -head.get_f64(120), // W_LONG
        phi: head.get_f64(72),   // S_LAT
    };

    let mut ur = Lp {
        lam: -head.get_f64(104), // E_LONG
        phi: head.get_f64(88),   // N_LAT
    };

    let mut del = Lp {
        lam: head.get_f64(152), // longitude interval
        phi: head.get_f64(136), // latitude interval
    };

    let lim = Lp {
        lam: (((ur.lam - ll.lam).abs() / del.lam + 0.5) + 1.).floor(),
        phi: (((ur.phi - ll.phi).abs() / del.phi + 0.5) + 1.).floor(),
    };

    // units are in seconds of degree.
    ll.lam *= SEC_TO_RAD;
    ll.phi *= SEC_TO_RAD;
    ur.lam *= SEC_TO_RAD;
    ur.phi *= SEC_TO_RAD;
    del.lam *= SEC_TO_RAD;
    del.phi *= SEC_TO_RAD;

    // Read matrix data
    let nrows = lim.phi as usize;
    let rowsize = lim.lam as usize;

    let gs_count = head.get_u32(168) as usize;
    if gs_count != nrows * rowsize {
        return Err(Error::InvalidNtv2GridFormat(ERR_GSCOUNT_NOT_MATCHING));
    }

    // Read grid data
    trace!(
        "Reading  data for grid {}:{}:{}",
        key,
        id.as_str(),
        lineage.as_str()
    );

    let mut buf = head.rebind::<16>();
    let mut cvs: Vec<Lp> = (0..gs_count)
        .map(|_| {
            buf.read(read)?;
            Ok(Lp {
                lam: SEC_TO_RAD * (buf.get_f32(0) as f64),
                phi: SEC_TO_RAD * (buf.get_f32(4) as f64),
            })
        })
        .collect::<Result<Vec<_>>>()?;

    // See https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NT111_V1_HARMEL_TransfoNTF-RGF93_FormatGrilleNTV2.pdf

    // In proj4, rows are stored in reverse order
    for i in 0..nrows {
        let offs = i * rowsize;
        cvs[offs..(offs + rowsize)].reverse();
    }

    let epsilon = (del.lam.abs() + del.phi.abs()) / 10_000.;

    catalog.add_grid(
        key.into(),
        Grid {
            id,
            lineage,
            ll,
            ur,
            del,
            lim,
            epsilon,
            cvs: cvs.into_boxed_slice(),
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nadgrids::Catalog;
    use crate::tests::setup;
    use std::env;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    macro_rules! fixture {
        ($name:expr) => {
            Path::new(&env::var("CARGO_MANIFEST_DIR").unwrap())
                .join("fixtures")
                .as_path()
                .join($name)
                .as_path()
        };
    }

    macro_rules! load_ntv2 {
        ($cat:expr, $name:expr) => {
            // Use a BufReader or efficiency
            let file = File::open(fixture!($name)).unwrap();
            let mut read = BufReader::new(file);
            read_ntv2($cat, $name, &mut read).unwrap();
        };
    }

    #[test]
    fn ntv2_100800401_gsb() {
        setup();

        let catalog = Catalog::default();
        load_ntv2!(&catalog, "100800401.gsb");

        let grids = catalog.find("100800401.gsb").unwrap().collect::<Vec<_>>();
        assert_eq!(grids.len(), 1);

        let grid = grids[0];
        assert!(grid.is_root());
        assert_eq!(grid.id.as_str(), "0INT2GRS");
        assert_eq!(grid.cvs.len(), 1591);
    }

    #[test]
    #[cfg(feature = "local_tests")]
    fn ntv2_bwta2017_gsb() {
        setup();

        let catalog = Catalog::default();
        load_ntv2!(&catalog, "BWTA2017.gsb");

        let grids = catalog.find("BWTA2017.gsb").unwrap().collect::<Vec<_>>();
        assert_eq!(grids.len(), 1);

        let grid = grids[0];
        assert!(grid.is_root());
        assert_eq!(grid.id.as_str(), "DHDN90  ");
        assert_eq!(grid.cvs.len(), 24514459);
    }
}
