//!
//! Handle Nadgrids
//!
use crate::errors::{Error, Result};

/// Nadgrid shift operator
pub trait NadgridShift: PartialEq {
    fn apply(&self, inverse: bool, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)>;
}

/// Nadgrid Provider trait
pub trait NadgridProvider {
    type ShiftItem: NadgridShift;

    fn find(&self, gridlist: &str) -> Result<Self::ShiftItem>;
}

//
// Implement a dummy grid shift operator and provider
//
use std::ops::ControlFlow;

#[derive(PartialEq)]
pub struct NullGridShift {}

impl NadgridShift for NullGridShift {
    fn apply(&self, _inverse: bool, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((x, y, z))
    }
}

/// Provider that return error if there is any mandatory
/// grid.
///
/// see https://proj.org/usage/transformation.html#the-null-grid
/// for the meaning of the null grid
pub struct NullGridProvider {}

impl NadgridProvider for NullGridProvider {
    type ShiftItem = NullGridShift;

    fn find(&self, gridlist: &str) -> Result<Self::ShiftItem> {
        // Parse the grid list and return an error
        // if there is any mandatory grid or the list is not terminated by
        // '@null'
        match gridlist.split(',').try_for_each(|s| {
            let s = s.trim();
            if s == "@null" || s == "null" {
                ControlFlow::Break(true)
            } else if s.starts_with('@') {
                // Optional grid
                ControlFlow::Continue(())
            } else {
                // Mand
                ControlFlow::Break(false)
            }
        }) {
            ControlFlow::Break(true) => Ok(NullGridShift {}),
            _ => Err(Error::NoNADGridAvailable),
        }
    }
}
