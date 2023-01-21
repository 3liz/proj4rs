//!
//! Overall coordinate system to coordinate system transformations
//! including reprojection and datum shifting
//!

use crate::consts::{EPS_12, FRAC_PI_2};
use crate::datum_transform::Datum;
use crate::errors::{Error, Result};
use crate::geocent::{geocentric_to_geodetic, geodetic_to_geocentric};
use crate::math::adjlon;
use crate::proj::{Axis, Proj};
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
/// Single point transform example:
///
/// ```rust
/// use proj4rs::transform::{transform, Transform};
/// use proj4rs::errors::Result;
///
/// pub struct Point {
///     x: f64,
///     y: f64,
///     z: f64,
/// }
///
/// impl Transform for Point {
///     fn transform_coordinates<F>(&mut self, mut f: F) -> Result<()>
///     where
///         F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>,
///     {
///         f(self.x, self.y, self.z).map(|(x, y, z)| {
///             self.x = x;
///             self.y = y;
///             self.z = z;
///         })
///     }
/// }
/// ```
///

pub trait Transform {
    fn transform_coordinates<F>(&mut self, f: F) -> Result<()>
    where
        F: FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>;
}

// ------------------
// Transformation
// ------------------

/// Select transformation direction
#[derive(PartialEq)]
pub enum Direction {
    Forward,
    Inverse,
}

use Direction::*;

/// The transformation function
pub fn transform<P>(src: &Proj, dst: &Proj, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
{
    if !src.has_inverse() {
        return Err(Error::NoInverseProjectionDefined);
    }

    if !dst.has_forward() {
        return Err(Error::NoForwardProjectionDefined);
    }

    adjust_axes(src, Inverse, points)?;
    geographic_to_cartesian(src, Inverse, points)?;
    projected_to_geographic(src, points)?;
    prime_meridian(src, Inverse, points)?;
    height_unit(src, Inverse, points)?;
    //geometric_to_orthometric(src, Inverse, points)?;

    datum_transform(src, dst, points)?;

    //geometric_to_orthometric(dst, Forward, points)?;
    height_unit(dst, Forward, points)?;
    prime_meridian(dst, Forward, points)?;
    geographic_to_cartesian(dst, Forward, points)?;
    geographic_to_projected(dst, points)?;
    //long_wrap(dst)?;
    adjust_axes(dst, Forward, points)?;

    Ok(())
}
// ---------------------------------
// Datum transformation
// ---------------------------------
fn datum_transform<P>(src: &Proj, dst: &Proj, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
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
// Projected to geographic
// ---------------------------------
fn projected_to_geographic<P>(p: &Proj, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
{
    // Nothing to do ?
    if p.is_latlong && !p.geoc
    /* && p.vto_meter == 1. */
    {
        return Ok(());
    }

    let (lam0, x0, y0) = (p.lam0(), p.x0(), p.y0());
    let (ra, one_es, to_meter) = (p.ellps.ra, p.ellps.one_es, p.to_meter);

    let pj_inv = p.inverse().ok_or(Error::NoInverseProjectionDefined)?;

    // Input points are cartesians
    // proj4 source: pj_inv.c
    points.transform_coordinates(|x, y, z| {
        if x.is_nan() || y.is_nan() {
            Err(Error::NanCoordinateValue)
        } else {
            // Inverse project
            let (mut lam, mut phi, z) = pj_inv(
                &p.projection,
                // descale and de-offset
                // z is not scaled since that
                // is handled by vto_meter before we get here
                (x * to_meter - x0) * ra,
                (y * to_meter - y0) * ra,
                z,
            )?;
            lam += lam0;
            if !p.over {
                lam = adjlon(lam);
            }
            if p.geoc && (phi.abs() - FRAC_PI_2).abs() > EPS_12 {
                phi = (one_es * phi.tan()).atan();
            }
            Ok((lam, phi, z))
        }
    })
}
// ---------------------------------
// Geographic to projected
// ---------------------------------
fn geographic_to_projected<P>(p: &Proj, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
{
    // Nothing to do ?
    //if (p.is_latlong && !p.geoc && p.vto_meter == 1.) || p.is_geocent {
    if (p.is_latlong && !p.geoc) || p.is_geocent {
        return Ok(());
    }

    let (lam0, x0, y0) = (p.lam0(), p.x0(), p.y0());
    let (a, rone_es, to_meter) = (p.ellps.a, p.ellps.rone_es, p.to_meter);

    let pj_fwd = p.forward().ok_or(Error::NoForwardProjectionDefined)?;
    let fr_meter = 1. / p.to_meter;

    // Input points are geographic
    // proj4 source: pj_fwd.c
    points.transform_coordinates(|lam, phi, z| {
        if lam.is_nan() || phi.is_nan() {
            Err(Error::NanCoordinateValue)
        } else {
            // Overrange check
            let t = phi.abs() - FRAC_PI_2;
            if t > EPS_12 || lam.abs() > 10. {
                Err(Error::CoordinateOutOfRange)
            } else {
                let (x, y, z) = pj_fwd(
                    &p.projection,
                    // ----
                    // lam
                    // ----
                    if !p.over {
                        adjlon(lam - lam0)
                    } else {
                        lam - lam0
                    },
                    // ---
                    // phi
                    // ---
                    if t.abs() <= EPS_12 {
                        if phi < 0. {
                            -FRAC_PI_2
                        } else {
                            FRAC_PI_2
                        }
                    } else if p.geoc {
                        (rone_es * phi.tan()).atan()
                    } else {
                        phi
                    },
                    // ---
                    // z
                    // ---
                    z,
                )?;
                // Rescale and offset
                Ok((fr_meter * (a * x + x0), fr_meter * (a * y + y0), z))
            }
        }
    })
}
// ---------------------------------
// Transform cartesian ("geocentric")
// source coordinates to lat/long,
// ---------------------------------
fn geographic_to_cartesian<P>(p: &Proj, dir: Direction, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
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
fn prime_meridian<P>(p: &Proj, dir: Direction, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
{
    let mut pm = p.from_greenwich;
    if pm == 0. || p.is_geocent || p.is_latlong {
        Ok(())
    } else {
        if dir == Forward {
            pm = -pm;
        }
        points.transform_coordinates(|x, y, z| Ok((x + pm, y, z)))
    }
}
// ---------------------
// Axis
// --------------------
fn adjust_axes<P>(p: &Proj, dir: Direction, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
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
fn normalize_axis<P: Transform + ?Sized>(axis: &Axis, points: &mut P) -> Result<()> {
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
fn denormalize_axis<P: Transform + ?Sized>(axis: &Axis, points: &mut P) -> Result<()> {
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
// ---------------------
// Adjust for vertical
// scale factor if needed
// --------------------
fn height_unit<P>(p: &Proj, dir: Direction, points: &mut P) -> Result<()>
where
    P: Transform + ?Sized,
{
    //if p.is_latlong {
    //    return Ok(());
    //}

    let fac = if dir == Forward {
        1. / p.vto_meter
    } else {
        p.vto_meter
    };

    if fac != 1.0 {
        points.transform_coordinates(|x, y, z| Ok((x, y, z * fac)))
    } else {
        Ok(())
    }
}
