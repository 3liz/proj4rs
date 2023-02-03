//!
//! Handle Nadgrids
//!
use crate::errors::{Error, Result};
use crate::transform::Direction;

mod grid;

#[cfg(feature = "multi-thread")]
mod catlg_mt;

#[cfg(feature = "multi-thread")]
use catlg_mt::{catalog, GridRef};

#[cfg(any(not(feature = "multi-thread"), target_arch = "wasm32"))]
mod catlg_st;

#[cfg(any(not(feature = "multi-thread"), target_arch = "wasm32"))]
use catlg_st::{catalog, GridRef};

use std::ops::ControlFlow;

/// NadGrids
///
/// Returned from the sequence
/// of nadgrids from projstring definition
#[derive(Debug)]
pub struct NadGrids(Vec<GridRef>);

impl PartialEq for NadGrids {
    fn eq(&self, other: &Self) -> bool {
        // Don't bother to compare all names
        self.0.is_empty() && other.0.is_empty()
    }
}

impl NadGrids {
    pub fn apply_shift(
        &self,
        dir: Direction,
        lam: f64,
        phi: f64,
        z: f64,
    ) -> Result<(f64, f64, f64)> {
        // Find the correct grid for an input
        match self.0.iter().find_map(|g| g.find_grid(lam, phi, z)) {
            Some(g) => g.nad_cvt(dir, lam, phi, z),
            None => Err(Error::PointOutsideNadShiftArea),
        }
    }

    /// Return a list of grids from the catalog
    pub fn new_grid_transform(names: &str) -> Result<Self> {
        // Parse the grid list and return an error
        // if there is any mandatory grid or the list is not terminated by
        // '@null'
        let mut v: Vec<GridRef> = vec![];

        match names.split(',').try_for_each(|s| {
            let s = s.trim();
            if s == "@null" || s == "null" {
                // Allow empty list
                // Mark also the end of parsing
                ControlFlow::Break(true)
            } else if let Some(s) = s.strip_prefix('@') {
                // Optional grid
                if let Some(g) = catalog::find_grid(s) {
                    v.push(g);
                }
                ControlFlow::Continue(())
            } else {
                // Mandatory grid
                if let Some(g) = catalog::find_grid(s) {
                    v.push(g);
                    ControlFlow::Continue(())
                } else {
                    ControlFlow::Break(false)
                }
            }
        }) {
            ControlFlow::Break(true) => Ok(Self(v)),
            ControlFlow::Break(false) => Err(Error::NadGridNotAvailable),
            _ => {
                if v.is_empty() {
                    Err(Error::NadGridNotAvailable)
                } else {
                    Ok(Self(v))
                }
            }
        }
    }
}
