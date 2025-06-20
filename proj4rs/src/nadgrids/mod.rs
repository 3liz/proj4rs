//!
//! Handle Nadgrids
//!
use crate::errors::{Error, Result};
use crate::transform::Direction;

mod catlg;
pub(crate) mod grid;

pub use catlg::{catalog, Catalog, GridRef};

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
mod header;

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
pub mod files;

use std::ops::ControlFlow;

pub use grid::Grid;

/// NadGrids
///
/// Returned from the sequence
/// of nadgrids from projstring definition
#[derive(Debug, Clone)]
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
        if self.0.is_empty() {
            return Ok((lam, phi, z));
        }

        // Find the correct (root)  grid for an input
        let mut iter = self.0.iter();
        let mut candidate = iter.find(|g| g.is_root() && g.matches(lam, phi, z));

        // Check for childs grid
        if let Some(grid) = candidate {
            let _ = iter.try_fold(grid, |grid, g| {
                if !g.is_child_of(grid) {
                    // No more childs, stop with the last candidate
                    ControlFlow::Break(())
                } else if g.matches(lam, phi, z) {
                    // Match, check for childs
                    candidate.replace(g);
                    ControlFlow::Continue(g)
                } else {
                    // Go next child
                    ControlFlow::Continue(grid)
                }
            });
        }

        match candidate {
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
                catalog::find_grids(s, &mut v);
                ControlFlow::Continue(())
            } else {
                // Mandatory grid
                if catalog::find_grids(s, &mut v) {
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

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}
