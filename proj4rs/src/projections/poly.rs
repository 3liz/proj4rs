//!
//! Implementation of the poly (American Polyconic) projection.
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
//! Parameters:
//!
//! * lat_0, lon_0, x_0, y_0
//!
//! Used by EPSG:5880 (SIRGAS 2000 / Brazil Polyconic), EPSG:5530, EPSG:29101 and EPSG:5472.
//!
//! See <https://proj.org/en/stable/operations/projections/poly.html>
//! and Snyder, *Map Projections: A Working Manual*, USGS PP 1395, pp. 124–137.
//!

use crate::errors::{Error, Result};
use crate::math::{Enfn, aasin, enfn, mlfn, msfn};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { poly }

const TOL: f64 = 1e-10;
const CONV: f64 = 1e-10;
const N_ITER: usize = 10;
const I_ITER: usize = 20;
const ITOL: f64 = 1e-12;

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    Ell {
        ml0: f64,
        es: f64,
        one_es: f64,
        en: Enfn,
    },
    Sph {
        phi0: f64,
    },
}

impl Projection {
    pub fn poly(p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        Ok(if p.ellps.is_ellipsoid() {
            let en = enfn(p.ellps.es);
            let (s0, c0) = p.phi0.sin_cos();
            Self::Ell {
                ml0: mlfn(p.phi0, s0, c0, en),
                es: p.ellps.es,
                one_es: p.ellps.one_es,
                en,
            }
        } else {
            Self::Sph { phi0: p.phi0 }
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Ell { ml0, es, en, .. } => {
                if phi.abs() <= TOL {
                    Ok((lam, -ml0, z))
                } else {
                    let (sp, cp) = phi.sin_cos();
                    let ms = if cp.abs() > TOL {
                        msfn(sp, cp, *es) / sp
                    } else {
                        0.
                    };
                    let e = lam * sp;
                    Ok((
                        ms * e.sin(),
                        (mlfn(phi, sp, cp, *en) - ml0) + ms * (1. - e.cos()),
                        z,
                    ))
                }
            }
            Self::Sph { phi0 } => {
                if phi.abs() <= TOL {
                    Ok((lam, -phi0, z))
                } else {
                    let cot = 1. / phi.tan();
                    let e = lam * phi.sin();
                    Ok((e.sin() * cot, phi - phi0 + cot * (1. - e.cos()), z))
                }
            }
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Ell {
                ml0,
                es,
                one_es,
                en,
            } => {
                let y = y + ml0;
                if y.abs() <= TOL {
                    return Ok((x, 0., z));
                }
                // Newton-Raphson on the meridional distance (Snyder, eq. 18-17 to 18-19).
                let r = y * y + x * x;
                let mut phi = y;
                let mut converged = false;
                for _ in 0..I_ITER {
                    let (sp, cp) = phi.sin_cos();
                    if cp.abs() < ITOL {
                        return Err(Error::CoordTransOutsideProjectionDomain);
                    }
                    let s2ph = sp * cp;
                    let w = (1. - es * sp * sp).sqrt();
                    let c = sp * w / cp;
                    let ml = mlfn(phi, sp, cp, *en);
                    let mlb = ml * ml + r;
                    let mlp = one_es / (w * w * w);
                    let dphi = (ml + ml + c * mlb - 2. * y * (c * ml + 1.))
                        / (es * s2ph * (mlb - 2. * y * ml) / c
                            + 2. * (y - ml) * (c * mlp - 1. / s2ph)
                            - mlp
                            - mlp);
                    phi += dphi;
                    if dphi.abs() <= ITOL {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                let sp = phi.sin();
                let lam = aasin(x * phi.tan() * (1. - es * sp * sp).sqrt())? / sp;
                Ok((lam, phi, z))
            }
            Self::Sph { phi0 } => {
                let y = phi0 + y;
                if y.abs() <= TOL {
                    return Ok((x, 0., z));
                }
                let b = x * x + y * y;
                let mut phi = y;
                let mut converged = false;
                for _ in 0..N_ITER {
                    let tp = phi.tan();
                    let dphi = (y * (phi * tp + 1.) - phi - 0.5 * (phi * phi + b) * tp)
                        / ((phi - y) / tp - 1.);
                    phi -= dphi;
                    if dphi.abs() <= CONV {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                let lam = aasin(x * phi.tan())? / phi.sin();
                Ok((lam, phi, z))
            }
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
    fn proj_poly_e() {
        let p = Proj::from_proj_string("+proj=poly +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (12.09, 47.73, 0.),
                (903237.329319470, 5359068.140851167, 0.),
            ),
            ((2.0, 1.0, 0.), (222605.285770237, 110642.194561440, 0.)),
            (
                (-3.0, -40.0, 0.),
                (-256133.208875131, -4433839.682716963, 0.),
            ),
            // On the equator
            ((10.0, 0.0, 0.), (1113194.907932736, 0.000000000, 0.)),
        ];

        // At 47.73° the 8th-degree meridional distance series (`mlfn`, "accurate to < 1e-5 m")
        // differs from PROJ 9 by about 1.1 µm, hence the tolerance.
        test_proj_forward(&p, &inputs, 1e-5);
        test_proj_inverse(&p, &inputs, 1e-5);
    }

    #[test]
    fn proj_poly_e_epsg_5880() {
        // EPSG:5880, SIRGAS 2000 / Brazil Polyconic
        let p = Proj::from_proj_string(
            "+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            // Brasília, Manaus, Porto Alegre, Rio Branco, Recife
            (
                (-47.9292, -15.7801, 0.),
                (5650400.250493755, 8245319.346777052, 0.),
            ),
            (
                (-60.0217, -3.119, 0.),
                (4330657.413116232, 9653201.624517998, 0.),
            ),
            (
                (-51.2177, -30.0346, 0.),
                (5268334.180255776, 6672789.93350167, 0.),
            ),
            (
                (-67.81, -9.9747, 0.),
                (3486203.259224749, 8865338.748656629, 0.),
            ),
            (
                (-34.877, -8.0476, 0.),
                (7107170.149865377, 9060846.324469948, 0.),
            ),
            // Origin, and the equator away from the central meridian
            (
                (-54.0, 0.0, 0.),
                (5000000.000000000, 10000000.000000000, 0.),
            ),
            ((-40.0, 0.0, 0.), (6558472.87110583, 10000000.000000000, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_poly_e_lat_0() {
        let p = Proj::from_proj_string("+proj=poly +lat_0=30 +lon_0=-96 +ellps=clrk66").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (-100.0, 40.0, 0.),
                (-341469.766841660, 1117048.629680726, 0.),
            ),
            ((-90.0, 25.0, 0.), (605513.446349582, -540635.893823085, 0.)),
            ((-96.0, 30.0, 0.), (0.000000000, 0.000000000, 0.)),
            (
                (-80.0, 0.0, 0.),
                (1781131.232828392, -3319933.299107925, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_poly_s() {
        let p = Proj::from_proj_string("+proj=poly +R_a +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (12.09, 47.73, 0.),
                (900068.906806543, 5374785.440031578, 0.),
            ),
            (
                (-3.0, -40.0, 0.),
                (-255349.948085921, -4449612.438589429, 0.),
            ),
            ((10.0, 0.0, 0.), (1111328.742043959, 0.000000000, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_poly_s_lat_0() {
        let p = Proj::from_proj_string("+proj=poly +lat_0=30 +lon_0=-96 +R=6371000").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (-100.0, 40.0, 0.),
                (-340606.678562030, 1119592.918486304, 0.),
            ),
            ((-90.0, 25.0, 0.), (604463.601577840, -542596.733702585, 0.)),
            (
                (-80.0, 0.0, 0.),
                (1779118.826312941, -3335847.799336762, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }
}
