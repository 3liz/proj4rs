//!
//! Transform adaptors
//!
use crate::errors::Result;
use crate::proj::Proj;
use crate::transform::{transform, Transform};

//
// Transform a 3-tuple
//
impl Transform for (f64, f64, f64) {
    fn transform_coordinates<F>(&mut self, f: &mut F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        (self.0, self.1, self.2) = f(self.0, self.1, self.2)?;
        Ok(())
    }
}

//
// Transform a 2-tuple
//
impl Transform for (f64, f64) {
    fn transform_coordinates<F>(&mut self, f: &mut F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        (self.0, self.1) = f(self.0, self.1, 0.).map(|(x, y, _)| (x, y))?;
        Ok(())
    }
}

/// Transform a 3-tuple
///
/// ```rust
/// use proj4rs::Proj;
/// use proj4rs::adaptors::transform_vertex_3d;
///
/// let dst = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();
/// let src = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();

/// let (x, y, z) = transform_vertex_3d(&src, &dst, (2.0, 1.0, 0.0)).unwrap();
/// ```
pub fn transform_vertex_3d(src: &Proj, dst: &Proj, pt: (f64, f64, f64)) -> Result<(f64, f64, f64)> {
    let mut pt_out = pt;
    transform(src, dst, &mut pt_out)?;
    Ok(pt_out)
}

/// Transform a 2-tuple
///
/// ```rust
/// use proj4rs::Proj;
/// use proj4rs::adaptors::transform_vertex_2d;
///
/// let dst = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();
/// let src = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();

/// let (x, y) = transform_vertex_2d(&src, &dst, (2.0, 1.0)).unwrap();
/// ```
#[inline(always)]
pub fn transform_vertex_2d(src: &Proj, dst: &Proj, pt: (f64, f64)) -> Result<(f64, f64)> {
    transform_vertex_3d(src, dst, (pt.0, pt.1, 0.)).map(|(x, y, _)| (x, y))
}

/// Transform x, y and z value
///
/// ```rust
/// use proj4rs::Proj;
/// use proj4rs::adaptors::transform_xyz;
///
/// let dst = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();
/// let src = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();

/// let (x, y, z) = transform_xyz(&src, &dst, 2.0, 1.0, 0.0).unwrap();
/// ```
#[inline(always)]
pub fn transform_xyz(src: &Proj, dst: &Proj, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
    transform_vertex_3d(src, dst, (x, y, z))
}

/// Transform x, y value
///
/// ```rust
/// use proj4rs::Proj;
/// use proj4rs::adaptors::transform_xy;
///
/// let dst = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();
/// let src = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();

/// let (x, y) = transform_xy(&src, &dst, 2.0, 1.0).unwrap();
/// ```
#[inline(always)]
pub fn transform_xy(src: &Proj, dst: &Proj, x: f64, y: f64) -> Result<(f64, f64)> {
    transform_xyz(src, dst, x, y, 0.).map(|(x, y, _)| (x, y))
}

//
// Transform an array of 3-tuple:
//
impl Transform for [(f64, f64, f64)] {
    fn transform_coordinates<F>(&mut self, f: &mut F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        self.iter_mut()
            .try_for_each(|xyz| xyz.transform_coordinates(f))
    }
}

//
// Transform an array of 2-tuple:
//
impl Transform for [(f64, f64)] {
    fn transform_coordinates<F>(&mut self, f: &mut F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
    {
        self.iter_mut()
            .try_for_each(|xy| xy.transform_coordinates(f))
    }
}
