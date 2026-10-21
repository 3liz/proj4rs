//!
//! Robinson
//!
//! ref: <https://proj.org/en/stable/operations/projections/robin.html>
//!
//! robin: "Robinson" "\n\tPCyl, Sph";
//!
//! Ported from PROJ `src/projections/robin.cpp`. Spherical only: on an
//! ellipsoid, the semi-major axis is used as the radius, as in PROJ.
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
use crate::errors::{Error, Result};
use crate::math::consts::{FRAC_PI_2, PI};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { robin }

// Coefficients are single precision, as in PROJ, and terms are based upon
// 5 degree intervals in degrees.
// See <https://trac.osgeo.org/proj/ticket/113>
type Coefs = (f32, f32, f32, f32);

#[rustfmt::skip]
const X: [Coefs; 19] = [
    (1.0, 2.2199e-17, -7.15515e-05, 3.1103e-06),
    (0.9986, -0.000482243, -2.4897e-05, -1.3309e-06),
    (0.9954, -0.00083103, -4.48605e-05, -9.86701e-07),
    (0.99, -0.00135364, -5.9661e-05, 3.6777e-06),
    (0.9822, -0.00167442, -4.49547e-06, -5.72411e-06),
    (0.973, -0.00214868, -9.03571e-05, 1.8736e-08),
    (0.96, -0.00305085, -9.00761e-05, 1.64917e-06),
    (0.9427, -0.00382792, -6.53386e-05, -2.6154e-06),
    (0.9216, -0.00467746, -0.00010457, 4.81243e-06),
    (0.8962, -0.00536223, -3.23831e-05, -5.43432e-06),
    (0.8679, -0.00609363, -0.000113898, 3.32484e-06),
    (0.835, -0.00698325, -6.40253e-05, 9.34959e-07),
    (0.7986, -0.00755338, -5.00009e-05, 9.35324e-07),
    (0.7597, -0.00798324, -3.5971e-05, -2.27626e-06),
    (0.7186, -0.00851367, -7.01149e-05, -8.6303e-06),
    (0.6732, -0.00986209, -0.000199569, 1.91974e-05),
    (0.6213, -0.010418, 8.83923e-05, 6.24051e-06),
    (0.5722, -0.00906601, 0.000182, 6.24051e-06),
    (0.5322, -0.00677797, 0.000275608, 6.24051e-06),
];

#[rustfmt::skip]
const Y: [Coefs; 19] = [
    (-5.20417e-18, 0.0124, 1.21431e-18, -8.45284e-11),
    (0.062, 0.0124, -1.26793e-09, 4.22642e-10),
    (0.124, 0.0124, 5.07171e-09, -1.60604e-09),
    (0.186, 0.0123999, -1.90189e-08, 6.00152e-09),
    (0.248, 0.0124002, 7.10039e-08, -2.24e-08),
    (0.31, 0.0123992, -2.64997e-07, 8.35986e-08),
    (0.372, 0.0124029, 9.88983e-07, -3.11994e-07),
    (0.434, 0.0123893, -3.69093e-06, -4.35621e-07),
    (0.4958, 0.0123198, -1.02252e-05, -3.45523e-07),
    (0.5571, 0.0121916, -1.54081e-05, -5.82288e-07),
    (0.6176, 0.0119938, -2.41424e-05, -5.25327e-07),
    (0.6769, 0.011713, -3.20223e-05, -5.16405e-07),
    (0.7346, 0.0113541, -3.97684e-05, -6.09052e-07),
    (0.7903, 0.0109107, -4.89042e-05, -1.04739e-06),
    (0.8435, 0.0103431, -6.4615e-05, -1.40374e-09),
    (0.8936, 0.00969686, -6.4636e-05, -8.547e-06),
    (0.9394, 0.00840947, -0.000192841, -4.2106e-06),
    (0.9761, 0.00616527, -0.000256, -4.2106e-06),
    (1.0, 0.00328947, -0.000319159, -4.2106e-06),
];

const FXC: f64 = 0.8487;
const FYC: f64 = 1.3523;
const C1: f64 = 11.459_155_902_616_464;
const RC1: f64 = 0.087_266_462_599_716_48;
const NODES: usize = 18;
const ONEEPS: f64 = 1.000001;
const EPS: f64 = 1e-10;
const MAX_ITER: usize = 100;

#[inline(always)]
fn v(c: &Coefs, z: f64) -> f64 {
    c.0 as f64 + z * (c.1 as f64 + z * (c.2 as f64 + z * c.3 as f64))
}

#[inline(always)]
fn dv(c: &Coefs, z: f64) -> f64 {
    c.1 as f64 + 2. * z * c.2 as f64 + z * z * 3. * c.3 as f64
}

#[derive(Debug, Clone)]
pub(crate) struct Projection {}

impl Projection {
    pub fn robin(_p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        Ok(Self {})
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if phi.is_nan() {
            return Err(Error::CoordinateOutOfRange);
        }
        let dphi = phi.abs();
        let i = ((dphi * C1 + 1e-15).floor() as usize).min(NODES);
        let dphi = (dphi - RC1 * i as f64).to_degrees();
        let x = v(&X[i], dphi) * FXC * lam;
        let y = v(&Y[i], dphi) * FYC;
        Ok((x, if phi < 0. { -y } else { y }, z))
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut lam = x / FXC;
        let mut phi = (y / FYC).abs();
        if phi >= 1. {
            // Simple pathologic cases
            if phi > ONEEPS {
                return Err(Error::CoordinateOutOfRange);
            }
            phi = FRAC_PI_2.copysign(y);
            lam /= X[NODES].0 as f64;
        } else {
            // In Y space, reduce to table interval
            if phi.is_nan() {
                return Err(Error::CoordinateOutOfRange);
            }
            let mut i = (phi * NODES as f64).floor() as usize;
            if i >= NODES {
                return Err(Error::CoordinateOutOfRange);
            }
            loop {
                if Y[i].0 as f64 > phi {
                    i -= 1;
                } else if Y[i + 1].0 as f64 <= phi {
                    i += 1;
                } else {
                    break;
                }
            }
            let t_ = &Y[i];
            // First guess, linear interpolation
            let mut t = 5. * (phi - t_.0 as f64) / (Y[i + 1].0 as f64 - t_.0 as f64);
            // Newton-Raphson
            let mut converged = false;
            for _ in 0..MAX_ITER {
                let t1 = (v(t_, t) - phi) / dv(t_, t);
                t -= t1;
                if t1.abs() < EPS {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err(Error::CoordinateOutOfRange);
            }
            phi = (5. * i as f64 + t).to_radians();
            if y < 0. {
                phi = -phi;
            }
            lam /= v(&X[i], t);
            if lam.abs() > PI {
                return Err(Error::CoordinateOutOfRange);
            }
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
    fn proj_robin_sphere() {
        let p = Proj::from_proj_string("+proj=robin +R=6370997").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((0., 0., 0.), (0., 0., 0.)),
            (
                (150., -30., 0.),
                (13589436.628231298, -3204965.794434718, 0.),
            ),
            (
                (-100., 60., 0.),
                (-7536475.401958603, 6328945.808890665, 0.),
            ),
            ((179., 89., 0.), (9109195.514557514, 8584445.796890883, 0.)),
            ((20., -90., 0.), (1004485.83992909, -8615499.2431, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs[2..], 1e-6);

        // The inverse is approximate near table nodes; same values as PROJ
        let inputs = [
            ((0., 0., 0.), (0., 0., 0.)),
            (
                (150.00000278883562, -30.00001449565707, 0.),
                (13589436.628231298, -3204965.794434718, 0.),
            ),
        ];
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_robin_ellipsoid() {
        // Spherical formulas applied with the major axis, as in PROJ
        let p = Proj::from_proj_string("+proj=robin +ellps=WGS84 +lon_0=150 +x_0=1000").unwrap();

        let inputs = [
            ((150., 0., 0.), (1000., 0., 0.)),
            ((-170., 45., 0.), (3387806.189311139, 4805073.646653474, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }
}
