//!
//! Implementations for some useful adaptors
//!
use crate::errors::Result;
use crate::proj::Proj;
use crate::transform::{transform, Transform};

//
// Transform a 3-tuple
//
impl Transform for (f64, f64, f64) {
    fn transform_coordinates<F>(&mut self, mut f: F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        (self.0, self.1, self.2) = f(self.0, self.1, self.2)?;
        Ok(())
    }
}

/// Transform a 3-tuple
pub fn transform_point_3d(src: &Proj, dst: &Proj, pt: (f64, f64, f64)) -> Result<(f64, f64, f64)> {
    let mut pt_out = pt;
    transform(src, dst, &mut pt_out)?;
    Ok(pt_out)
}

/// Transform a 2-tuple
#[inline(always)]
pub fn transform_point_2d(src: &Proj, dst: &Proj, pt: (f64, f64)) -> Result<(f64, f64)> {
    transform_point_3d(src, dst, (pt.0, pt.1, 0.)).map(|(x, y, _)| (x, y))
}

/// Transform x, y and z value
#[inline(always)]
pub fn transform_xyz(src: &Proj, dst: &Proj, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
    transform_point_3d(src, dst, (x, y, z))
}

/// Transform x, y value (z is set to 0)
#[inline(always)]
pub fn transform_xy(src: &Proj, dst: &Proj, x: f64, y: f64) -> Result<(f64, f64)> {
    transform_xyz(src, dst, x, y, 0.).map(|(x, y, _)| (x, y))
}

//
// Transform an array of 3-tuple:
//
impl Transform for [(f64, f64, f64)] {
    fn transform_coordinates<F>(&mut self, mut f: F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        self.iter_mut().try_for_each(|(x, y, z)| {
            (*x, *y, *z) = f(*x, *y, *z)?;
            Ok(())
        })
    }
}

/// Transform an array of 3 tuple
#[inline(always)]
pub fn transform_point_array(src: &Proj, dst: &Proj, pts: &mut [(f64, f64, f64)]) -> Result<()> {
    transform(src, dst, pts)
}
