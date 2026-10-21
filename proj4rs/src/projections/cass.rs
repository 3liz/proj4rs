//!
//! Cassini
//!
//! ref: <https://proj.org/en/stable/operations/projections/cass.html>
//!
//! cass: "Cassini" "\n\tCyl, Sph&Ell";
//!
//! Ported from PROJ `src/projections/cass.cpp`, with `+hyperbolic`. The
//! ellipsoidal inverse is refined by Newton iteration on the forward, as
//! PROJ does with `pj_generic_inverse_2d`.
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
use crate::math::{
    Enfn,
    consts::{FRAC_PI_2, PI},
    enfn, inv_mlfn, mlfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { cass }

const C1: f64 = 1. / 6.;
const C2: f64 = 1. / 120.;
const C3: f64 = 1. / 24.;
const C4: f64 = 1. / 3.;
const C5: f64 = 1. / 15.;

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    phi0: f64,
    es: f64,
    en: Enfn,
    m0: f64,
    hyperbolic: bool,
}

#[allow(non_snake_case)]
impl Projection {
    pub fn cass(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let es = p.ellps.es;
        let en = enfn(es);
        let m0 = if p.ellps.is_ellipsoid() {
            mlfn(p.phi0, p.phi0.sin(), p.phi0.cos(), en)
        } else {
            0.
        };
        Ok(Self {
            phi0: p.phi0,
            es,
            en,
            m0,
            hyperbolic: params.check_option("hyperbolic")?,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (x, y) = if self.es == 0. {
            self.s_forward(lam, phi)
        } else {
            self.e_forward(lam, phi)
        };
        Ok((x, y, z))
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (lam, phi) = if self.es == 0. {
            self.s_inverse(x, y)
        } else {
            self.e_inverse(x, y)?
        };
        Ok((lam, phi, z))
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }

    //
    // Sphere
    //

    fn s_forward(&self, lam: f64, phi: f64) -> (f64, f64) {
        (
            (phi.cos() * lam.sin()).asin(),
            phi.tan().atan2(lam.cos()) - self.phi0,
        )
    }

    fn s_inverse(&self, x: f64, y: f64) -> (f64, f64) {
        let dd = y + self.phi0;
        ((x.tan()).atan2(dd.cos()), (dd.sin() * x.cos()).asin())
    }

    //
    // Ellipsoid
    //

    fn e_forward(&self, lam: f64, phi: f64) -> (f64, f64) {
        let es = self.es;
        let (sinphi, cosphi) = phi.sin_cos();
        let M = mlfn(phi, sinphi, cosphi, self.en);

        let nu_square = 1. / (1. - es * sinphi * sinphi);
        let nu = nu_square.sqrt();
        let tanphi = phi.tan();
        let T = tanphi * tanphi;
        let A = lam * cosphi;
        let C = es * (cosphi * cosphi) / (1. - es);
        let A2 = A * A;

        let x = nu * A * (1. - A2 * T * (C1 + (8. - T + 8. * C) * A2 * C2));
        let mut y = M - self.m0 + nu * tanphi * A2 * (0.5 + (5. - T + 6. * C) * A2 * C3);
        if self.hyperbolic {
            let rho = nu_square * (1. - es) * nu;
            y -= y * y * y / (6. * rho * nu);
        }
        (x, y)
    }

    fn e_inverse(&self, x: f64, y: f64) -> Result<(f64, f64)> {
        let es = self.es;
        let phi1 = inv_mlfn(self.m0 + y, self.en);
        let tanphi1 = phi1.tan();
        let T1 = tanphi1 * tanphi1;
        let sinphi1 = phi1.sin();
        let nu1_square = 1. / (1. - es * sinphi1 * sinphi1);
        let nu1 = nu1_square.sqrt();
        let rho1 = nu1_square * (1. - es) * nu1;
        let D = x / nu1;
        let D2 = D * D;
        let phi = phi1 - (nu1 * tanphi1 / rho1) * D2 * (0.5 - (1. + 3. * T1) * D2 * C3);
        let lam = D * (1. + T1 * D2 * (-C4 + (1. + 3. * T1) * D2 * C5)) / phi1.cos();

        self.e_inverse_refine(x, y, lam, phi, 1e-12)
    }

    // Newton iteration on the forward, from PROJ `pj_generic_inverse_2d`
    fn e_inverse_refine(
        &self,
        x: f64,
        y: f64,
        mut lam: f64,
        mut phi: f64,
        tol: f64,
    ) -> Result<(f64, f64)> {
        let (mut dlam_x, mut dlam_y, mut dphi_x, mut dphi_y) = (0., 0., 0., 0.);
        for i in 0..15 {
            let (xa, ya) = self.e_forward(lam, phi);
            let (dx, dy) = (xa - x, ya - y);
            if dx.abs() < tol && dy.abs() < tol {
                return Ok((lam, phi));
            }

            if i == 0 || dx.abs() > 1e-6 || dy.abs() > 1e-6 {
                // Jacobian, only while not close to the result
                let dl = if lam > 0. { -1e-6 } else { 1e-6 };
                let (x2, y2) = self.e_forward(lam + dl, phi);
                let (dx_lam, dy_lam) = ((x2 - xa) / dl, (y2 - ya) / dl);

                let dp = if phi > 0. { -1e-6 } else { 1e-6 };
                let (x2, y2) = self.e_forward(lam, phi + dp);
                let (dx_phi, dy_phi) = ((x2 - xa) / dp, (y2 - ya) / dp);

                let det = dx_lam * dy_phi - dx_phi * dy_lam;
                if det != 0. {
                    dlam_x = dy_phi / det;
                    dlam_y = -dx_phi / det;
                    dphi_x = -dy_lam / det;
                    dphi_y = dx_lam / det;
                }
            }

            // Limit the correction to avoid overshoots from a bad first guess
            lam -= (dx * dlam_x + dy * dlam_y).clamp(-0.3, 0.3);
            lam = lam.clamp(-PI, PI);
            phi -= (dx * dphi_x + dy * dphi_y).clamp(-0.3, 0.3);
            phi = phi.clamp(-FRAC_PI_2, FRAC_PI_2);
        }
        Err(Error::CoordinateOutOfRange)
    }
}

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_cass_sphere() {
        let p = Proj::from_proj_string("+proj=cass +R=6370997 +lat_0=-42 +lon_0=147").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (288859.731012250, 1330553.868361174, 0.)),
            ((140., -50., 0.), (-499590.934807034, -913014.620776347, 0.)),
            ((100., 10., 0.), (-5123251.540627883, 6282074.446714829, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_cass_ellipsoid() {
        let p = Proj::from_proj_string(
            "+proj=cass +ellps=WGS84 +lat_0=-42 +lon_0=147 +x_0=500000 +y_0=10000000",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (500000., 10000000., 0.)),
            ((150., -30., 0.), (789425.749226268, 11327731.424806831, 0.)),
            ((140., -50., 0.), (-1136.129708443, 9087261.414570767, 0.)),
            (
                (155., -20., 0.),
                (1336856.319812124, 12419140.321031865, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_cass_hyperbolic() {
        let p = Proj::from_proj_string(concat!(
            "+proj=cass +lat_0=-16.25 +lon_0=179.33333333333 +x_0=50640 +y_0=67295 ",
            "+a=6378306.3696 +b=6356571.996 +hyperbolic"
        ))
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((179.33333333333, -16.25, 0.), (50640., 67295., 0.)),
            ((179., -16.5, 0.), (15050.862850907, 39601.982048825, 0.)),
            ((-179.8, -17., 0.), (142930.300404808, -15899.851477167, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-9);
    }
}
