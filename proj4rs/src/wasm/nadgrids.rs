//!
//! WASM provider for nadgrids
//!
//! Use JS Dataview for passing nadgrids definition
//!
use crate::errors::Error;
use crate::math::consts::SEC_TO_RAD;
use crate::nadgrids::grid::{Grid, GridId, Lp};
use crate::nadgrids::catalog;
use js_sys::DataView;
use wasm_bindgen::prelude::*;

const ERR_INVALID_HEADER: &str = "Wrong header";
const ERR_GSCOUNT_NOT_MATCHING: &str = "GS COUNT not matching";

const HEADER_SIZE: usize = 11 * 16;

/// Read a binary NTv2 from Dataview.
///
/// Note: only NTv2 file format are supported.
#[wasm_bindgen]
pub fn add_nadgrid(key: &str, view: &DataView) -> Result<(), JsError> {
    // Check endianess
    let is_le = view.get_int32_endian(8, true) == 11;

    // Read NTv2 overview header
    let nfields = view.get_int32_endian(8, is_le);
    if nfields != 11 {
        return Err(Error::InvalidNtv2GridFormat(ERR_INVALID_HEADER).into());
    }

    let nsubgrids = view.get_int32_endian(40, is_le) as usize;

    // Read subsequent grids
    (0..nsubgrids).try_fold(HEADER_SIZE, |offset, _i| {
        read_subgrid(view, offset, is_le).and_then(|grid| {
            let offs = offset + grid.gs_count() * 16 + HEADER_SIZE;
            catalog::add_grid(key.into(), grid)?;
            Ok(offs)
        })
    })?;
    Ok(())
}

fn read_subgrid(view: &DataView, offset: usize, is_le: bool) -> Result<Grid, Error> {
    match view
        .buffer()
        .slice_with_end(offset as u32, offset as u32 + 8)
        .as_string()
    {
        Some(s) if s == "SUB_NAME" => Ok(()),
        _ => Err(Error::InvalidNtv2GridFormat(ERR_INVALID_HEADER)),
    }?;

    // SUB_NAME
    let id = GridId::from((
        view.get_uint32_endian(offset + 4, is_le),
        view.get_uint32_endian(offset + 8, is_le),
    ));

    // PARENT
    let mut lineage = GridId::from((
        view.get_uint32_endian(offset + 24, is_le),
        view.get_uint32_endian(offset + 24 + 4, is_le),
    ));

    if lineage.as_str() == "NONE" {
        lineage = GridId::root();
    }

    let mut ll = Lp {
        lam: -view.get_float64_endian(120 + offset, is_le), // W_LONG
        phi: view.get_float64_endian(72 + offset, is_le),   // S_LAT
    };

    let ur = Lp {
        lam: -view.get_float64_endian(104 + offset, is_le), // E_LONG
        phi: view.get_float64_endian(88 + offset, is_le),   // N_LAT
    };

    let mut del = Lp {
        lam: view.get_float64_endian(152 + offset, is_le), // longitude interval
        phi: view.get_float64_endian(136 + offset, is_le), // latitude interval
    };

    let lim = Lp {
        lam: (((ur.lam - ll.lam).abs() / del.lam + 0.5) + 1.).floor(),
        phi: (((ur.phi - ll.phi).abs() / del.phi + 0.5) + 1.).floor(),
    };

    // units are in seconds of degree.
    ll.lam *= SEC_TO_RAD;
    ll.phi *= SEC_TO_RAD;
    del.lam *= SEC_TO_RAD;
    del.phi *= SEC_TO_RAD;

    // Read matrix data
    let nrows = lim.phi as usize;
    let rowsize = lim.lam as usize;

    let gs_count = view.get_int32_endian(168 + offset, is_le) as usize;
    if gs_count != nrows * rowsize {
        return Err(Error::InvalidNtv2GridFormat(ERR_GSCOUNT_NOT_MATCHING));
    }

    let cvsoffset = offset + HEADER_SIZE;
    let mut cvs: Vec<Lp> = (0..gs_count)
        .map(|i| Lp {
            lam: SEC_TO_RAD * (view.get_float32_endian(cvsoffset + i * 16, is_le) as f64),
            phi: SEC_TO_RAD * (view.get_float32_endian(cvsoffset + i * 16 + 4, is_le) as f64),
        })
        .collect();

    // See https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NT111_V1_HARMEL_TransfoNTF-RGF93_FormatGrilleNTV2.pdf

    // In proj4, rows are stored in reverse order
    for i in 0..nrows {
        let offs = i * rowsize;
        cvs[offs..(offs + rowsize)].reverse();
    }

    let epsilon = (del.lam.abs() + del.phi.abs()) / 10_000.;

    Ok(Grid {
        id,
        lineage,
        ll,
        ur,
        del,
        lim,
        epsilon,
        cvs: cvs.into_boxed_slice(),
    })
}
