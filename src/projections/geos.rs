//!
//! Geostationary Satellite View (Azimuthal)
//!
//! proj: geos
//!
//! ----------------------------------
//! Required
//! +h=<value> 
//! Height of the view point above the Earth and must be in the same units
//! as the radius of the sphere or semimajor axis of the ellipsoid.
//!
//! Optional
//! +sweep=<axis>¶
//! Sweep angle axis of the viewing instrument. Valid options are "x" and "y".
//! Defaults to "y".
//! 
//! lon_0: the reference longitude
//! ----------------------------------


// Projection stub
super::projection! { geos }

use crate::errors::{Error, Result};
use crate::math::consts::FRAC_PI_2;
use crate::parameters::ParamList;
use crate::proj::ProjData;

#[derive(Debug, Clone)]
pub(crate) struct Ell {
    radius_p: f64,
    radius_p2: f64,
    radius_p_inv2: f64,
    radius_g: f64,
    radius_g_1: f64,
    c: f64,
    flip_axis: i32,
}

#[derive(Debug, Clone)]
pub(crate) struct Sph {
    radius_g: f64,
    radius_g_1: f64,
    c: f64,
    flip_axis: i32,
}

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    Ell(Ell),
    Sph(Sph),
}

use Projection::*;

impl Projection {
    pub fn geos(p: &mut ProjData, params: &ParamList) -> Result<Self> {
      
        let h:f64 = params
            .try_value("h")?
            .ok_or_else(|| Error::InputStringError("Parameter h is required."))?;

        let flip_axis = 0; // TODO catch "sweep=x" or "sweep=y" (default)

        let radius_g_1 = h / p.ellps.a;
        if radius_g_1 <= 0. || radius_g_1 > 1e10 {
            return Err(Error::InvalidParameterValue("Invalid value for h."));
        }

        let radius_g = 1. + radius_g_1;
        let c = radius_g * radius_g - 1.0;

        if p.ellps.is_ellipsoid() {
            Ok(Ell(Ell {
                radius_p: (p.ellps.one_es).sqrt(),
                radius_p2: p.ellps.one_es,
                radius_p_inv2: p.ellps.rone_es,
                radius_g,
                radius_g_1,
                c,
                flip_axis,
            }))
        } else {
            Ok(Sph(Sph {
                radius_g,
                radius_g_1,
                c,
                flip_axis,
            }))
        }
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Ell(e) => e.forward(lam, phi, z),
            Sph(s) => s.forward(lam, phi, z),
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Ell(e) => e.inverse(x, y, z),
            Sph(s) => s.inverse(x, y, z),
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

// ---------------
// Ellipsoidal
// ---------------
#[rustfmt::skip]
impl Ell {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {

        // Calculation of geocentric latitude.
        let phie = (self.radius_p2 * phi.tan()).atan();

        // Calculation of the three components of the vector from satellite to
        // position on earth surface (long,lat).
        let cos_phi = phie.cos(); 
        let sin_phi = phie.sin(); 
        let r = self.radius_p / (self.radius_p * cos_phi).hypot(sin_phi);
    
        let vx = r * lam.cos() * cos_phi;
        let vy = r * lam.sin() * cos_phi;
        let vz = r * sin_phi;

        // Check visibility.
        if ((self.radius_g - vx) * vx - vy * vy - vz * vz * self.radius_p_inv2) < 0. {
            return Err(Error::CoordTransOutsideProjectionDomain);
        }

        // Calculation based on view angles from satellite.
        let tmp = self.radius_g - vx;
        if self.flip_axis == 1 {
            Ok((self.radius_g_1 * (vy / vz.hypot(tmp)).atan(),
                self.radius_g_1 * (vz / tmp).atan(),
                z))
        } else {
            Ok((self.radius_g_1 * (vy / tmp).atan(),
                self.radius_g_1 * (vz / vy.hypot(tmp)).atan(),
                z))
        }

    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Setting three components of vector from satellite to position.
        let mut vx = -1.0;
        let mut vy;
        let mut vz;

        if self.flip_axis == 1 {
            vz = (y / self.radius_g_1).tan();
            vy = (x / self.radius_g_1).tan() * 1.0_f64.hypot(vz);
        } else {
            vy = (x / self.radius_g_1).tan();
            vz = (y / self.radius_g_1).tan() * 1.0_f64.hypot(vy);
        }
    
        // Calculation of terms in cubic equation and determinant.
        let mut aa = vz / self.radius_p;
        aa = vy * vy + aa* aa + vx * vx;
        let b = 2.0 * self.radius_g * vx;
        let det = b * b - 4.0 * aa * self.c;

        if det < 0. {
            return Err(Error::CoordTransOutsideProjectionDomain);   
        }
    
        // Calculation of three components of vector from satellite to position.
        let k = (-b - det.sqrt()) / (2. * aa);
        vx = self.radius_g + k * vx;
        vy *= k;
        vz *= k;
    
        // Calculation of longitude and latitude.
        let lam = vy.atan2(vx);
        Ok((lam, 
            (self.radius_p_inv2 * ((vz * lam.cos() / vx).atan())),// .tan()).atan(), 
            z))
    }
}

// ---------------
// Spherical
// ---------------
impl Sph {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Fail if our longitude is more than 90 degrees from the
        // central meridian since the results are essentially garbage.
        if !(-FRAC_PI_2..=FRAC_PI_2).contains(&lam) {
            return Err(Error::LatOrLongExceedLimit);
        }

        let cosphi = phi.cos();
        let vx = cosphi * lam.cos();
        let vy = cosphi * lam.sin();
        let vz = phi.sin();
        let tmp = self.radius_g - vx;

        if self.flip_axis == 1 {
            Ok((
                self.radius_g_1 * (vy / vz.hypot(tmp)).atan(),
                self.radius_g_1 * (vz / tmp).atan(),
                z,
            ))
        } else {
            Ok((
                self.radius_g_1 * (vy / tmp).atan(),
                self.radius_g_1 * (vz / vy.hypot(tmp)).atan(),
                z,
            ))
        }
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Setting three components of vector from satellite to position.
        let mut vx = -1.0;
        let mut vy;
        let mut vz;
        if self.flip_axis == 1 {
            vz = (y / self.radius_g_1).tan();
            vy = (x / self.radius_g_1).tan() * (1.0 + vz * vz).sqrt();
        } else {
            vy = (x / self.radius_g_1).tan();
            vz = (y / self.radius_g_1).tan() * (1.0 + vy * vy).sqrt();
        }

        // Calculation of terms in cubic equation and determinant.
        let a = vy * vy + vz * vz + vx * vx;
        let b = 2.0 * self.radius_g * vx;

        let det = b * b - 4.0 * a * self.c;
        if det < 0. {
            return Err(Error::CoordTransOutsideProjectionDomain);
        }

        // Calculation of three components of vector from satellite to position.
        let k = (-b - det.sqrt()) / (2. * a);
        vx = self.radius_g + k * vx;
        vy *= k;
        vz *= k;

        // Calculation of longitude and latitude.
        let lam = vy.atan2(vx);

        Ok((lam, (vz * lam.cos() / vx).atan(), z))
    }
}

#[cfg(test)]
mod tests {
    use crate::math::consts::{EPS_10, EPS_7};
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_geos_el() {
        let p = Proj::from_proj_string("+proj=geos +lon_0=0 +h=35785782.858 +x_0=0 +y_0=0 +a=6378160 +b=6356775 +units=m +no_defs")
            .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            //
            (
                (18.763481601401576, 9.204293875870595, 0.),
                (2000000.0, 1000000.0, 0.),
            ),
            (
                (18.763481601401576, -9.204293875870595, 0.),
                (2000000.0, -1000000.0, 0.),
            ),
            (
                (-18.763481601401576, 9.204293875870595, 0.),
                (-2000000.0, 1000000.0, 0.),
            ),
            (
                (-18.763481601401576, -9.204293875870595, 0.),
                (-2000000.0, -1000000.0, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1.0e-2);
    }

    #[test]
    fn proj_geos_sp() {
        let p = Proj::from_proj_string("+proj=geos +lon_0=0 +h=35785833.8833").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (18.763554109081273, 9.204326881322723, 0.),
                (2000000.0, 1000000.0, 0.),
            ),
            (
                (18.763554109081273, -9.204326881322723, 0.),
                (2000000.0, -1000000.0, 0.),
            ),
            (
                (-18.763554109081273, 9.204326881322723, 0.),
                (-2000000.0, 1000000.0, 0.),
            ),
            (
                (-18.763554109081273, -9.204326881322723, 0.),
                (-2000000.0, -1000000.0, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1.0e-8);
        test_proj_inverse(&p, &inputs, 1.0e-2);
    }
}
