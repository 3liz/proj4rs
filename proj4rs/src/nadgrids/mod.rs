//!
//! Handle Nadgrids
//!
use crate::errors::{Error, Result};
use crate::transform::Direction;

mod catlg;
mod header;

pub mod formats;

pub(crate) mod grid;

pub use catlg::{Catalog, GridRef, catalog};

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
pub mod files;

use std::ops::ControlFlow;

pub use grid::Grid;



/// NadGrids
///
/// Returned from the sequence
/// of nadgrids from projstring definition
#[derive(Debug, Clone)]
pub struct NadGrids {
    grids: Vec<GridRef>,
    nullable: bool,
}


impl PartialEq for NadGrids {
    fn eq(&self, other: &Self) -> bool {
        // Compare references
        self.nullable == other.nullable
        && self.grids.len() == other.grids.len()
        && self.grids.iter().zip(&other.grids).all(|(g1, g2)| std::ptr::eq(*g1, *g2))
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
        if self.grids.is_empty() {
            return Ok((lam, phi, z));
        }

        // Find the correct (root)  grid for an input
        let mut iter = self.grids.iter();
        let mut candidate = iter.find(|g| g.is_root() && g.matches(lam, phi, z));

        // Check for childs grid
        if let Some(grid) = candidate {
            let _ = iter.fold(grid, |grid, g| {
                if !g.is_child_of(grid) {
                    // Skip it
                    grid
                } else if g.matches(lam, phi, z) {
                    // Match, check for childs
                    candidate.replace(g);
                    g
                } else {
                    // Go next child
                    grid
                }
            });
        }

        match candidate {
            Some(g) => g.nad_cvt(dir, lam, phi, z),
            None if self.nullable => Ok((lam, phi, z)),
            None => Err(Error::PointOutsideNadShiftArea),
        }
    }

    /// Return a list of grids from the catalog
    pub fn new_grid_transform(names: &str) -> Result<Self> {
        // Parse the grid list and return an error
        // if there is any mandatory grid or the list is not terminated by
        // '@null'
        let mut v: Vec<GridRef> = vec![];
        let mut nullable = false;

        match names.split(',').try_for_each(|s| {
            let s = s.trim();
            if s == "@null" || s == "null" {
                // Identity grid
                // Mark also the end of parsing
                nullable = true;
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
            ControlFlow::Break(true) => Ok(Self { grids: v, nullable }),
            ControlFlow::Break(false) => Err(Error::NadGridNotAvailable),
            _ => {
                if v.is_empty() {
                    Err(Error::NadGridNotAvailable)
                } else {
                    Ok(Self { grids: v, nullable })
                }
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.grids.is_empty()
    }
}
