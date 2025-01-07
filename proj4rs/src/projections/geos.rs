//!
//! From proj/geos.cpp
//!
//! See also <https://proj.org/operations/projections/geos.html>
//!
//! Geostationnary Satellite View
//!
//! The geos projection pictures how a geostationary satellite scans the earth at regular scanning angle intervals.
//!
//! Copyright (c) 2004   Gerald I. Evenden
//! Copyright (c) 2012   Martin Raspaud
//!
//!  See also (section 4.4.3.2):
//!    https://www.cgms-info.org/documents/pdf_cgms_03.pdf
//!
//! Permission is hereby granted, free of charge, to any person obtaining
//! a copy of this software and associated documentation files (the
//! "Software"), to deal in the Software without restriction, including
//! without limitation the rights to use, copy, modify, merge, publish,
//! distribute, sublicense, and/or sell copies of the Software, and to
//! permit persons to whom the Software is furnished to do so, subject to
//! the following conditions:
//!
//! The above copyright notice and this permission notice shall be
//! included in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//!
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { geos }

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    El(Ell),
    Sp(Sph),
}

impl Projection {
    pub fn geos(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let h: f64 = params
            .try_value("h")?
            .ok_or(Error::InputStringError("Missing parameter 'h'"))?;
        let flip_axis: bool = params
            .try_value::<&str>("sweep")
            .and_then(|sweep| match sweep {
                None => Ok(false),
                Some("y") => Ok(false),
                Some("x") => Ok(true),
                Some(_other) => Err(Error::InvalidParameterValue(
                    "sweep require only 'x' or 'y' value",
                )),
            })?;

        let radius_g_1 = h / p.ellps.a;
        if radius_g_1 <= 0. || radius_g_1 >= 1.0e10 {
            return Err(Error::InvalidParameterValue("Invalid value for 'h'."));
        }

        let radius_g = 1. + radius_g_1;
        let c = radius_g * radius_g - 1.0;

        if p.ellps.is_ellipsoid() {
            Ok(Self::El(Ell {
                radius_p: p.ellps.one_es.sqrt(),
                radius_p2: p.ellps.one_es,
                radius_p_inv2: p.ellps.rone_es,
                radius_g,
                radius_g_1,
                c,
                flip_axis,
            }))
        } else {
            Ok(Self::Sp(Sph {
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
            Self::El(p) => p.forward(lam, phi, z),
            Self::Sp(p) => p.forward(lam, phi, z),
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::El(p) => p.inverse(x, y, z),
            Self::Sp(p) => p.inverse(x, y, z),
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

//
// Ellipsoid
//

#[derive(Debug, Clone)]
pub(crate) struct Ell {
    radius_p: f64,
    radius_p2: f64,
    radius_p_inv2: f64,
    radius_g: f64,
    radius_g_1: f64,
    c: f64,
    flip_axis: bool,
}

impl Ell {
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Calculation of geocentric latitude.
        // g_phi = (self.radius_p2 * phi.tan()).atan();

        // Calculation of the three components of the vector from satellite to
        // position on earth surface (long,lat).
        let (sin_phi, cos_phi) = (self.radius_p2 * phi.tan()).atan().sin_cos();
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
        if self.flip_axis {
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

        if self.flip_axis {
            vz = (y / self.radius_g_1).tan();
            vy = (x / self.radius_g_1).tan() * 1.0_f64.hypot(vz);
        } else {
            vy = (x / self.radius_g_1).tan();
            vz = (y / self.radius_g_1).tan() * 1.0_f64.hypot(vy);
        }

        // Calculation of terms in cubic equation and determinant.
        let mut a = vz / self.radius_p;
        a = vy * vy + a * a + vx * vx;
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
        Ok((
            lam,
            (self.radius_p_inv2 * ((vz * lam.cos() / vx).atan())),
            z,
        ))
    }
}

//
// Spherical
//

#[derive(Debug, Clone)]
pub(crate) struct Sph {
    radius_g: f64,
    radius_g_1: f64,
    c: f64,
    flip_axis: bool,
}

impl Sph {
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut tmp = phi.cos();
        let vx = tmp * lam.cos();
        let vy = tmp * lam.sin();
        let vz = phi.sin();

        tmp = self.radius_g - vx;

        if self.flip_axis {
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
        if self.flip_axis {
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
