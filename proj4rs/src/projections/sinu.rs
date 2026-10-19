//!
//! Implementation of the sinu (Sinusoidal, Sanson-Flamsteed) projection.
//! Original author (Proj library) Author:   Gerald Evenden
//!
//!
//! Copyright (c) 1995, Gerald Evenden
//!
//! Permission is hereby granted, free of charge, to any person obtaining a
//! copy of this software and associated documentation files (the "Software"),
//! to deal in the Software without restriction, including without limitation
//! the rights to use, copy, modify, merge, publish, distribute, sublicense,
//! and/or sell copies of the Software, and to permit persons to whom the
//! Software is furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included
//! in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//! DEALINGS IN THE SOFTWARE.
//!
//! Ported from PROJ `src/projections/gn_sinu.cpp`, its `sinu` case (`m = 0`, `n = 1`), with
//! its ellipsoidal and spherical forms.
//!
//! Parameters:
//!
//! * lon_0, x_0, y_0
//!
//! Used by MODIS land products (sphere, R = 6371007.181), ESRI:54008 (World_Sinusoidal),
//! ESRI:53008 (Sphere_Sinusoidal) and ESRI:102011 (Africa_Sinusoidal).
//!
//! See <https://proj.org/en/stable/operations/projections/sinu.html>
//! and Snyder, *Map Projections: A Working Manual*, USGS PP 1395, pp. 243–248.
//!

use crate::errors::{Error, Result};
use crate::math::consts::{EPS_10, FRAC_PI_2};
use crate::math::{Enfn, enfn, inv_mlfn, mlfn};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { sinu }

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    Ell { es: f64, en: Enfn },
    Sph,
}

impl Projection {
    pub fn sinu(p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        Ok(if p.ellps.is_ellipsoid() {
            Self::Ell {
                es: p.ellps.es,
                en: enfn(p.ellps.es),
            }
        } else {
            Self::Sph
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Ell { es, en } => {
                let (s, c) = phi.sin_cos();
                Ok((lam * c / (1. - es * s * s).sqrt(), mlfn(phi, s, c, *en), z))
            }
            Self::Sph => Ok((lam * phi.cos(), phi, z)),
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (lam, phi) = match self {
            Self::Ell { es, en } => {
                let phi = inv_mlfn(y, *es, *en)?;
                let s = phi.abs();
                if s < FRAC_PI_2 {
                    let s = phi.sin();
                    (x * (1. - es * s * s).sqrt() / phi.cos(), phi)
                } else if (s - EPS_10) < FRAC_PI_2 {
                    (0., phi)
                } else {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
            }
            Self::Sph => {
                if y.abs() > FRAC_PI_2 + EPS_10 {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                let c = y.cos();
                (if c.abs() > EPS_10 { x / c } else { 0. }, y)
            }
        };
        if lam.is_finite() {
            Ok((lam, phi, z))
        } else {
            Err(Error::CoordTransOutsideProjectionDomain)
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

//============
// Tests
//============

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.5.1 (through pyproj 3.7.1).

    #[test]
    fn proj_sinu_e() {
        let p = Proj::from_proj_string("+proj=sinu +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (12.09, 47.73, 0.),
                (906918.099039101, 5288406.917509883, 0.),
            ),
            ((2.0, 1.0, 0.), (222605.299539466, 110574.388554153, 0.)),
            (
                (-3.0, -40.0, 0.),
                (-256181.570877594, -4429529.030236588, 0.),
            ),
            // On the equator
            ((10.0, 0.0, 0.), (1113194.907932736, 0.000000000, 0.)),
        ];

        // As for `poly`: the 8th-degree meridional distance series differs from PROJ 9 by about
        // a micrometre at mid latitudes.
        test_proj_forward(&p, &inputs, 1e-5);
        test_proj_inverse(&p, &inputs, 1e-5);
    }

    #[test]
    fn proj_sinu_e_esri_102011() {
        // ESRI:102011, Africa_Sinusoidal: a central meridian on WGS 84
        let p = Proj::from_proj_string(
            "+proj=sinu +lon_0=15 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            // Nairobi, Oslo, Sorriso (Brazil), and the central meridian on the equator
            (
                (36.82, -1.29, 0.),
                (2428379.790083543, -142641.057809881, 0.),
            ),
            (
                (10.75, 59.91, 0.),
                (-237793.838437986, 6644045.782491288, 0.),
            ),
            (
                (-55.72, -12.55, 0.),
                (-7685629.2129637, -1387927.939345077, 0.),
            ),
            ((15.0, 0.0, 0.), (0.000000000, 0.000000000, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-5);
        test_proj_inverse(&p, &inputs, 1e-5);
    }

    #[test]
    fn proj_sinu_e_false_origin() {
        let p =
            Proj::from_proj_string("+proj=sinu +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS80")
                .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            // Brasília, Manaus, and the origin
            (
                (-47.9292, -15.7801, 0.),
                (5650490.259000941, 8254690.364816656, 0.),
            ),
            (
                (-60.0217, -3.119, 0.),
                (4330653.76517124, 9655115.414839646, 0.),
            ),
            ((-54.0, 0.0, 0.), (5000000.0, 10000000.0, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-5);
        test_proj_inverse(&p, &inputs, 1e-5);
    }

    #[test]
    fn proj_sinu_s_modis() {
        // MODIS's grid: a sphere of radius 6371007.181 m
        let p = Proj::from_proj_string("+proj=sinu +R=6371007.181").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (-55.72, -12.55, 0.),
                (-6047750.815866172, -1395497.902306987, 0.),
            ),
            (
                (36.82, -1.29, 0.),
                (4093164.15381663, -143441.617049881, 0.),
            ),
            ((10.75, 59.91, 0.), (599298.755817259, 6661695.56392124, 0.)),
            ((179.0, 0.0, 0.), (19903914.303820767, 0.000000000, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_sinu_s_modis_geokeys() {
        // The same sphere as a MODIS GeoTIFF's GeoKeys give it (semi-major = semi-minor)
        let p = Proj::from_proj_string(
            "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
        )
        .unwrap();

        let inputs = [(
            (-55.72, -12.55, 0.),
            (-6047750.815866172, -1395497.902306987, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }
}
