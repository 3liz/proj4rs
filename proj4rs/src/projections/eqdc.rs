//!
//! Equidistant Conic
//!
//! ref: <https://proj.org/en/stable/operations/projections/eqdc.html>
//!
//! eqdc: "Equidistant Conic" "\n\tConic, Sph&Ell\n\tlat_1= lat_2=";
//!
//! Ported from PROJ `src/projections/eqdc.cpp`.
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
    consts::{EPS_10, FRAC_PI_2},
    enfn, inv_mlfn, mlfn, msfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { eqdc }

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    n: f64,
    rho0: f64,
    c: f64,
    en: Enfn,
    ellips: bool,
}

impl Projection {
    pub fn eqdc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let phi1 = params.try_angular_value("lat_1")?.unwrap_or(0.);
        let phi2 = params.try_angular_value("lat_2")?.unwrap_or(0.);

        if phi1.abs() > FRAC_PI_2 {
            return Err(Error::InvalidParameterValue("|lat_1| should be <= 90°"));
        }
        if phi2.abs() > FRAC_PI_2 {
            return Err(Error::InvalidParameterValue("|lat_2| should be <= 90°"));
        }
        if (phi1 + phi2).abs() < EPS_10 {
            return Err(Error::InvalidParameterValue(
                "|lat_1 + lat_2| should be > 0",
            ));
        }

        let es = p.ellps.es;
        let en = enfn(es);
        let (mut sinphi, mut cosphi) = phi1.sin_cos();
        let mut n = sinphi;
        let secant = (phi1 - phi2).abs() >= EPS_10;
        let ellips = p.ellps.is_ellipsoid();

        let (c, rho0);
        if ellips {
            let m1 = msfn(sinphi, cosphi, es);
            let ml1 = mlfn(phi1, sinphi, cosphi, en);
            if secant {
                (sinphi, cosphi) = phi2.sin_cos();
                let ml2 = mlfn(phi2, sinphi, cosphi, en);
                if ml1 == ml2 {
                    return Err(Error::InvalidParameterValue("Eccentricity too close to 1"));
                }
                n = (m1 - msfn(sinphi, cosphi, es)) / (ml2 - ml1);
                if n == 0. {
                    return Err(Error::InvalidParameterValue(
                        "Invalid value for eccentricity",
                    ));
                }
            }
            c = ml1 + m1 / n;
            rho0 = c - mlfn(p.phi0, p.phi0.sin(), p.phi0.cos(), en);
        } else {
            if secant {
                n = (cosphi - phi2.cos()) / (phi2 - phi1);
            }
            if n == 0. {
                return Err(Error::InvalidParameterValue("lat_1 + lat_2 should be > 0"));
            }
            c = phi1 + phi1.cos() / n;
            rho0 = c - p.phi0;
        }

        Ok(Self {
            n,
            rho0,
            c,
            en,
            ellips,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rho = self.c
            - if self.ellips {
                mlfn(phi, phi.sin(), phi.cos(), self.en)
            } else {
                phi
            };
        let (sinl, cosl) = (lam * self.n).sin_cos();
        Ok((rho * sinl, self.rho0 - rho * cosl, z))
    }

    pub fn inverse(&self, mut x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut y = self.rho0 - y;
        let mut rho = x.hypot(y);
        if rho != 0. {
            if self.n < 0. {
                rho = -rho;
                x = -x;
                y = -y;
            }
            let mut phi = self.c - rho;
            if self.ellips {
                phi = inv_mlfn(phi, self.en);
            }
            Ok((x.atan2(y) / self.n, phi, z))
        } else {
            Ok((0., FRAC_PI_2.copysign(self.n), z))
        }
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
    use crate::errors::Error;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_eqdc_sphere() {
        let p = Proj::from_proj_string("+proj=eqdc +R=6370997 +lat_1=-20 +lat_2=-40 +lon_0=135")
            .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((135., 0., 0.), (0., 0., 0.)),
            (
                (150., -30., 0.),
                (1418501.215579857, -3428347.067879791, 0.),
            ),
            (
                (110., -10., 0.),
                (-2831229.798074735, -1420434.57547921, 0.),
            ),
            (
                (170., -60., 0.),
                (2270179.643113833, -7019307.634750781, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_eqdc_ellipsoid() {
        let p = Proj::from_proj_string(
            "+proj=eqdc +ellps=GRS80 +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((-96., 40., 0.), (0., 0., 0.)),
            (
                (-120., 30., 0.),
                (-2193479.005798156, -818253.437336989, 0.),
            ),
            ((-70., 55., 0.), (1588170.056599415, 1896314.962503969, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_eqdc_tangent() {
        let p = Proj::from_proj_string("+proj=eqdc +ellps=WGS84 +lat_1=50 +lat_2=50 +lon_0=10")
            .unwrap();

        let inputs = [
            ((10., 50., 0.), (0., 5540847.041684148, 0.)),
            ((30., 60., 0.), (1122746.182793978, 6805084.847435765, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_eqdc_invalid() {
        assert!(matches!(
            Proj::from_proj_string("+proj=eqdc +R=6370997 +lat_1=20 +lat_2=-20"),
            Err(Error::InvalidParameterValue(_))
        ));
        assert!(matches!(
            Proj::from_proj_string("+proj=eqdc +R=6370997 +lat_1=91 +lat_2=20"),
            Err(Error::InvalidParameterValue(_))
        ));
    }
}
