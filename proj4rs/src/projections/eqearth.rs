//!
//! Equal Earth
//!
//! ref: <https://proj.org/en/stable/operations/projections/eqearth.html>
//!
//! eqearth: "Equal Earth" "\n\tPCyl, Sph&Ell";
//!
//! Ported from PROJ `src/projections/eqearth.cpp`.
//!
//! Copyright (c) 2018 Bojan Savric
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
use crate::errors::{Error, Result};
use crate::math::Authalic;
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { eqearth }

// Polynomial coefficients
const A1: f64 = 1.340264;
const A2: f64 = -0.081106;
const A3: f64 = 0.000893;
const A4: f64 = 0.003796;

// 90° latitude on a sphere with radius 1
const MAX_Y: f64 = 1.3173627591574;
const EPS: f64 = 1e-11;
const MAX_ITER: usize = 12;

#[inline(always)]
fn m() -> f64 {
    3f64.sqrt() / 2.
}

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    rqda: f64,
    /// Ellipsoidal case only
    auth: Option<Authalic>,
}

impl Projection {
    pub fn eqearth(p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        let el = &p.ellps;
        Ok(if el.is_ellipsoid() {
            let auth = Authalic::new(el);
            Self {
                // Authalic radius divided by major axis
                rqda: (0.5 * auth.qp()).sqrt(),
                auth: Some(auth),
            }
        } else {
            Self {
                rqda: 1.,
                auth: None,
            }
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut sbeta = phi.sin();
        if let Some(auth) = &self.auth {
            // Sine of the authalic latitude
            sbeta = (auth.q(sbeta) / auth.qp()).clamp(-1., 1.);
        }

        let psi = (m() * sbeta).asin();
        let psi2 = psi * psi;
        let psi6 = psi2 * psi2 * psi2;

        let x = lam * psi.cos() / (m() * (A1 + 3. * A2 * psi2 + psi6 * (7. * A3 + 9. * A4 * psi2)));
        let y = psi * (A1 + A2 * psi2 + psi6 * (A3 + A4 * psi2));
        Ok((x * self.rqda, y * self.rqda, z))
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let x = x / self.rqda;
        let y = (y / self.rqda).clamp(-MAX_Y, MAX_Y);

        // Newton-Raphson
        let mut yc = y;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let y2 = yc * yc;
            let y6 = y2 * y2 * y2;
            let f = yc * (A1 + A2 * y2 + y6 * (A3 + A4 * y2)) - y;
            let fder = A1 + 3. * A2 * y2 + y6 * (7. * A3 + 9. * A4 * y2);
            let tol = f / fder;
            yc -= tol;
            if tol.abs() < EPS {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(Error::CoordinateOutOfRange);
        }

        let y2 = yc * yc;
        let y6 = y2 * y2 * y2;
        let lam = m() * x * (A1 + 3. * A2 * y2 + y6 * (7. * A3 + 9. * A4 * y2)) / yc.cos();
        let mut phi = (yc.sin() / m()).asin();
        if let Some(auth) = &self.auth {
            phi = auth.lat_inverse(phi);
        }
        Ok((lam, phi, z))
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_eqearth_sphere() {
        let p = Proj::from_proj_string("+proj=eqearth +R=6370997 +lon_0=150").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((150., 0., 0.), (0., 0., 0.)),
            ((-170., 45., 0.), (3284190.904098491, 5480529.66516747, 0.)),
            (
                (100., -80., 0.),
                (-2955925.250371363, -8207189.666485644, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_eqearth_ellipsoid() {
        let p = Proj::from_proj_string("+proj=eqearth +ellps=WGS84").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((0., 0., 0.), (0., 0., 0.)),
            (
                (150., -30., 0.),
                (13448009.907913297, -3764325.426892128, 0.),
            ),
            (
                (-100., 60., 0.),
                (-7225128.360577337, 6924035.137552897, 0.),
            ),
            ((179., 89., 0.), (10164195.230873333, 8391009.94247709, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-10);
    }
}
