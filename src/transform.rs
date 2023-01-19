//!
//! Overall coordinate system to coordinate system transformations
//! including reprojection and datum shifting
//!

use crate::datum_transform::Datum;
use crate::errors::Result;
use crate::geocent::{geocentric_to_geodetic, geodetic_to_geocentric};
use crate::nadgrids::NadgridShift;
use crate::proj::{Axis, Projection};

///
/// Transform trait
///
/// This allow transform to be agnostic in the coordinate's
/// implementation details (useful for collections of coordinates)
///
/// The closure will return error if processing of the coordinate fail.
/// If The method may return an error, then the whole processing
/// stop. Indeed, the strategy to stop or continue on error
/// is left to the `ApplyTransform` implementation.
///
///
/// Example:
///
/// ```rust
///
/// use proj4js::transform::{transform, , Transform};
/// use crate::errors::Result;
///
/// pub struct Point {
///     x: f64,
///     y: f64,
///     z: f64,
/// }
///
/// impl Transform for Point {
///     fn transform_coordinates<F>(&mut self, mut f: F) -> Result<()>
///        where
///            F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
///        {
///            f(self.x, self.y, self.z).map(|(x, y, z)| {
///                self.x = x;
///                self.y = y;
///                self.z = z;
///            })
///        }
///    }
/// }
/// ```
///
pub trait Transform {
    fn transform_coordinates<F>(&mut self, f: F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>;
}

/// Select transformation direction
#[derive(PartialEq)]
pub enum Direction {
    Forward,
    Inverse,
}

use Direction::*;

// ------------------
// Projection
// ------------------
pub fn transform<N, P>(src: &Projection<N>, dst: &Projection<N>, points: &mut P) -> Result<()>
where
    N: NadgridShift,
    P: Transform,
{
    adjust_axes(src, Inverse, points)?;
    geographic_to_cartesian(src, Inverse, points)?;
    //projected_to_geographic(src, points)?;
    prime_meridian(src, Inverse, points)?;
    //height_unit(src, Inverse, points)?;
    //geometric_to_orthometric(src, Inverse, points)?;

    datum_transform(src, dst, points)?;

    //geometric_to_orthometric(dst, Forward, points)?;
    //height_unit(dst, Forward, points)?;
    prime_meridian(dst, Forward, points)?;
    geographic_to_cartesian(dst, Forward, points)?;
    //geographic_to_projected(dst, Forward, points)?;
    //long_wrap(dst)?;
    adjust_axes(dst, Forward, points)?;

    Ok(())
}
// ---------------------------------
// Datum transformation
// ---------------------------------
fn datum_transform<N, P>(src: &Projection<N>, dst: &Projection<N>, points: &mut P) -> Result<()>
where
    N: NadgridShift,
    P: Transform,
{
    let src_datum = &src.datum;
    let dst_datum = &dst.datum;

    // Return true if the datums are identical is respect
    // to datum transformation.
    // As of PROJ 4.6.0 behavior, we prevent datum transformation
    // if either the source or destination are of an unknown datum type.
    if src_datum.no_datum() || dst_datum.no_datum() || src_datum.is_identical_to(dst_datum) {
        return Ok(());
    }

    points.transform_coordinates(|x, y, z| Datum::transform(src_datum, dst_datum, x, y, z))
}
// ---------------------------------
// Transform cartesian ("geocentric")
// source coordinates to lat/long,
// ---------------------------------
fn geographic_to_cartesian<N, P>(p: &Projection<N>, dir: Direction, points: &mut P) -> Result<()>
where
    N: NadgridShift,
    P: Transform,
{
    // Nothing to do
    if !p.is_geocent {
        return Ok(());
    }

    let (a, b, es) = (p.datum.a, p.datum.b, p.datum.es);
    let fac = p.to_meter;

    if fac != 1.0 {
        match dir {
            Forward => points.transform_coordinates(|x, y, z| {
                geodetic_to_geocentric(x, y, z, a, es).map(|(x, y, z)| (x * fac, y * fac, z * fac))
            }),
            Inverse => points.transform_coordinates(|x, y, z| {
                geocentric_to_geodetic(x * fac, y * fac, z * fac, a, es, b)
            }),
        }
    } else {
        match dir {
            Forward => {
                points.transform_coordinates(|x, y, z| geodetic_to_geocentric(x, y, z, a, es))
            }
            Inverse => {
                points.transform_coordinates(|x, y, z| geocentric_to_geodetic(x, y, z, a, es, b))
            }
        }
    }
}
// --------------------------
// Prime meridian adjustement
// -------------------------
fn prime_meridian<N, P>(p: &Projection<N>, dir: Direction, points: &mut P) -> Result<()>
where
    N: NadgridShift,
    P: Transform,
{
    let mut pm = p.pm;
    if pm == 0.  || p.is_geocent || p.is_latlong {
        return Ok(());
    } 

    if dir == Forward {
        pm = -pm;
    }

    points.transform_coordinates(|x, y, z| Ok((x + pm, y, z))) 
}
// ---------------------
// Axis
// --------------------
fn adjust_axes<N, P>(p: &Projection<N>, dir: Direction, points: &mut P) -> Result<()>
where
    N: NadgridShift,
    P: Transform,
{
    if !p.normalized_axis() {
        match dir {
            Forward => denormalize_axis(&p.axis, points),
            Inverse => normalize_axis(&p.axis, points),
        }
    } else {
        Ok(())
    }
}

// Normalize axis
fn normalize_axis<P: Transform>(axis: &Axis, points: &mut P) -> Result<()> {
    points.transform_coordinates(|x, y, z| {
        let (mut x_out, mut y_out, mut z_out) = (x, y, z);
        axis.iter().enumerate().for_each(|(i, axe)| {
            let value = match i {
                1 => x,
                2 => y,
                _ => z,
            };
            match axe {
                b'e' => x_out = value,
                b'w' => x_out = -value,
                b'n' => y_out = value,
                b's' => y_out = -value,
                b'u' => z_out = value,
                b'd' => z_out = -value,
                // This is unreachable because the way
                // axis have been initialized (see the `proj` module)
                _ => unreachable!(),
            }
        });
        Ok((x_out, y_out, z_out))
    })
}

// Denormalize axis
fn denormalize_axis<P: Transform>(axis: &Axis, points: &mut P) -> Result<()> {
    points.transform_coordinates(|x, y, z| {
        let (mut x_out, mut y_out, mut z_out) = (x, y, z);
        axis.iter().enumerate().for_each(|(i, axe)| {
            let value = match axe {
                b'e' => x,
                b'w' => -x,
                b'n' => y,
                b's' => -y,
                b'u' => z,
                b'd' => -z,
                // See above
                _ => unreachable!(),
            };
            match i {
                1 => x_out = value,
                2 => y_out = value,
                _ => z_out = value,
            }
        });
        Ok((x_out, y_out, z_out))
    })
}
