//!
//! Gnomonic
//!
//! ref: <https://proj.org/en/stable/operations/projections/gnom.html>
//!
//! gnom: "Gnomonic" "\n\tAzi, Sph&Ell";
//!
//! Ported from PROJ `src/projections/gnom.cpp`. The ellipsoidal form is
//! Karney's (2013) geodesic gnomonic.
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
use proj4rs_geodesic::{Flags, Geodesic, Mask};

use crate::errors::{Error, Result};
use crate::math::consts::{EPS_10, FRAC_PI_2};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { gnom }

#[allow(non_camel_case_types)]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq)]
enum Mode {
    N_POLE,
    S_POLE,
    EQUIT,
    OBLIQ,
}

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    phi0: f64,
    sinph0: f64,
    cosph0: f64,
    mode: Mode,
    // Ellipsoidal form when set
    g: Option<Box<Geodesic>>,
}

impl Projection {
    pub fn gnom(p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        let phi0 = p.phi0;
        let (sinph0, cosph0) = phi0.sin_cos();

        let mode = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
            if phi0 < 0. {
                Mode::S_POLE
            } else {
                Mode::N_POLE
            }
        } else if phi0.abs() < EPS_10 {
            Mode::EQUIT
        } else {
            Mode::OBLIQ
        };

        let g = if p.ellps.is_sphere() {
            None
        } else {
            Some(Box::new(Geodesic::new(1., p.ellps.f)))
        };
        Ok(Self {
            phi0,
            sinph0,
            cosph0,
            mode,
            g,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if let Some(g) = &self.g {
            return self.e_forward(g, lam, phi, z);
        }
        self.s_forward(lam, phi, z)
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if let Some(g) = &self.g {
            return self.e_inverse(g, x, y, z);
        }
        self.s_inverse(x, y, z)
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

    fn s_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let coslam = lam.cos();

        let cosz = match self.mode {
            Mode::EQUIT => cosphi * coslam,
            Mode::OBLIQ => self.sinph0 * sinphi + self.cosph0 * cosphi * coslam,
            Mode::S_POLE => -sinphi,
            Mode::N_POLE => sinphi,
        };

        if cosz <= EPS_10 {
            return Err(Error::CoordinateOutOfRange);
        }

        let r = 1. / cosz;
        let x = r * cosphi * lam.sin();
        let y = r * match self.mode {
            Mode::EQUIT => sinphi,
            Mode::OBLIQ => self.cosph0 * sinphi - self.sinph0 * cosphi * coslam,
            Mode::N_POLE => -cosphi * coslam,
            Mode::S_POLE => cosphi * coslam,
        };
        Ok((x, y, z))
    }

    fn s_inverse(&self, mut x: f64, mut y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rh = x.hypot(y);
        let mut phi = rh.atan();
        let sinz = phi.sin();
        let cosz = (1. - sinz * sinz).sqrt();

        if rh.abs() <= EPS_10 {
            return Ok((0., self.phi0, z));
        }

        match self.mode {
            Mode::OBLIQ => {
                phi = cosz * self.sinph0 + y * sinz * self.cosph0 / rh;
                phi = if phi.abs() >= 1. {
                    FRAC_PI_2.copysign(phi)
                } else {
                    phi.asin()
                };
                y = (cosz - self.sinph0 * phi.sin()) * rh;
                x *= sinz * self.cosph0;
            }
            Mode::EQUIT => {
                phi = y * sinz / rh;
                phi = if phi.abs() >= 1. {
                    FRAC_PI_2.copysign(phi)
                } else {
                    phi.asin()
                };
                y = cosz * rh;
                x *= sinz;
            }
            Mode::S_POLE => phi -= FRAC_PI_2,
            Mode::N_POLE => {
                phi = FRAC_PI_2 - phi;
                y = -y;
            }
        }
        Ok((x.atan2(y), phi, z))
    }

    //
    // Ellipsoid (Karney, 2013, "Algorithms for geodesics", section 8)
    //

    fn e_forward(&self, g: &Geodesic, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let r = g.gen_inverse(
            self.phi0.to_degrees(),
            0.,
            phi.to_degrees(),
            lam.to_degrees(),
            Mask::AZIMUTH | Mask::REDUCED_LENGTH | Mask::GEODESIC_SCALE,
        );
        if r.scale12 <= 0. {
            return Err(Error::CoordinateOutOfRange);
        }
        let rho = r.m12 / r.scale12;
        let (sinazi, cosazi) = r.azi1.to_radians().sin_cos();
        Ok((rho * sinazi, rho * cosazi, z))
    }

    fn e_inverse(&self, g: &Geodesic, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        const NUMIT: usize = 10;
        let eps = 0.01 * f64::EPSILON.sqrt();

        // Azimuth clockwise from north
        let azi0 = x.atan2(y).to_degrees();
        let mut rho = x.hypot(y);
        let mut s = rho.atan();
        let little = rho <= 1.;
        if !little {
            rho = 1. / rho;
        }
        let l = g.line(
            self.phi0.to_degrees(),
            0.,
            azi0,
            Mask::LATITUDE
                | Mask::LONGITUDE
                | Mask::DISTANCE_IN
                | Mask::REDUCED_LENGTH
                | Mask::GEODESIC_SCALE,
        );

        let mut trip = false;
        let (mut lat1, mut lon1) = (0., 0.);
        for _ in 0..NUMIT {
            let p = l.gen_position(
                Flags::NONE,
                s,
                Mask::LATITUDE | Mask::LONGITUDE | Mask::REDUCED_LENGTH | Mask::GEODESIC_SCALE,
            );
            (lat1, lon1) = (p.lat2, p.lon2);
            let (m, mm) = (p.m12, p.scale12);
            if trip {
                break;
            }
            // If little, solve rho(s) = rho with drho(s)/ds = 1/M^2
            // else solve 1/rho(s) = 1/rho with d(1/rho(s))/ds = -1/m^2
            let ds = if little {
                (m - rho * mm) * mm
            } else {
                (rho * m - mm) * m
            };
            s -= ds;
            // Also stop on NaN
            if ds.abs() < eps || ds.is_nan() {
                trip = true;
            }
        }
        if trip {
            Ok((lon1.to_radians(), lat1.to_radians(), z))
        } else {
            Err(Error::CoordinateOutOfRange)
        }
    }
}

//============
// Tests
//============

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_gnom_s_obliq() {
        let p = Proj::from_proj_string("+proj=gnom +R=6370997 +lat_0=-42 +lon_0=147").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (295478.296841106, 1350242.102235714, 0.)),
            (
                (100., -60., 0.),
                (-2797139.686113018, -3177553.796137307, 0.),
            ),
            ((-170., -30., 0.), (4672935.288165815, 413276.766432783, 0.)),
            ((147., 10., 0.), (0., 8154504.298877209, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_s_equit() {
        let p = Proj::from_proj_string("+proj=gnom +R=6370997 +lat_0=0 +lon_0=0").unwrap();

        let inputs = [
            ((10., 20., 0.), (1123378.665112688, 2354625.320024798, 0.)),
            ((-80., 5., 0.), (-36131719.458938956, 3209881.153643807, 0.)),
            ((45., -45., 0.), (6370996.999999998, -9009950.363238299, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_s_poles() {
        let p = Proj::from_proj_string("+proj=gnom +R=6370997 +lat_0=90").unwrap();
        let inputs = [
            ((10., 60., 0.), (638729.541960903, -3622415.238975614, 0.)),
            (
                (-120., 10., 0.),
                (-31290986.933853678, 18065859.729469474, 0.),
            ),
            ((170., 85., 0.), (96789.759982481, 548922.006113809, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string("+proj=gnom +R=6370997 +lat_0=-90 +lon_0=30").unwrap();
        let inputs = [
            ((10., -60., 0.), (-1258051.610002067, 3456468.390989212, 0.)),
            (
                (-120., -10., 0.),
                (-18065859.729469493, -31290986.933853671, 0.),
            ),
            ((170., -85., 0.), (358283.394027932, -426985.521999635, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_e_obliq() {
        let p = Proj::from_proj_string("+proj=gnom +ellps=WGS84 +lat_0=-42 +lon_0=147").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (296026.117484902, 1347275.543766303, 0.)),
            (
                (100., -60., 0.),
                (-2807859.531655353, -3182754.956885839, 0.),
            ),
            ((-170., -30., 0.), (4683812.890027311, 407828.140511455, 0.)),
            ((147., 10., 0.), (0., 8094012.723442596, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_e_equit() {
        let p = Proj::from_proj_string("+proj=gnom +ellps=WGS84 +lat_0=0 +lon_0=0").unwrap();

        let inputs = [
            ((10., 20., 0.), (1124713.736254139, 2342041.514079806, 0.)),
            ((-80., 5., 0.), (-37065392.748346493, 3279116.307989170, 0.)),
            ((45., -45., 0.), (6388314.046392297, -8985889.731100047, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_e_poles() {
        let p = Proj::from_proj_string("+proj=gnom +ellps=WGS84 +lat_0=90").unwrap();
        let inputs = [
            ((10., 60., 0.), (640943.892175542, -3634973.443090167, 0.)),
            (
                (-120., 10., 0.),
                (-30616287.588145059, 17676321.880602539, 0.),
            ),
            ((170., 85., 0.), (97221.705750164, 551371.692293126, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string("+proj=gnom +ellps=WGS84 +lat_0=-90 +lon_0=30").unwrap();
        let inputs = [
            ((10., -60., 0.), (-1262413.028520589, 3468451.289885645, 0.)),
            (
                (-120., -10., 0.),
                (-17676321.880602557, -30616287.588145047, 0.),
            ),
            ((170., -85., 0.), (359882.313125461, -428891.039578382, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_e_k0() {
        let p = Proj::from_proj_string(
            "+proj=gnom +ellps=GRS80 +lat_0=40 +lon_0=-100 +k_0=0.9996 +x_0=1000 +y_0=-500",
        )
        .unwrap();

        let inputs = [
            ((-100., 40., 0.), (1000., -500., 0.)),
            ((-90., 45., 0.), (795033.992842076, 605635.247739543, 0.)),
            (
                (-120., 20., 0.),
                (-2286446.231335942, -2164835.129677983, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_gnom_hidden_hemisphere() {
        use crate::transform::transform;

        for (src, dst) in [
            (
                "+proj=longlat +R=6370997",
                "+proj=gnom +R=6370997 +lat_0=-42 +lon_0=147",
            ),
            (
                "+proj=longlat +R=6370997",
                "+proj=gnom +R=6370997 +lat_0=90",
            ),
            (
                "+proj=longlat +ellps=WGS84",
                "+proj=gnom +ellps=WGS84 +lat_0=-42 +lon_0=147",
            ),
            (
                "+proj=longlat +ellps=WGS84",
                "+proj=gnom +ellps=WGS84 +lat_0=90",
            ),
        ] {
            let src = Proj::from_proj_string(src).unwrap();
            let p = Proj::from_proj_string(dst).unwrap();
            // Antipode of the centre, or a point on the far hemisphere
            let mut pt = ((-33f64).to_radians(), (42f64).to_radians(), 0.);
            if dst.contains("lat_0=90") {
                pt.1 = (-10f64).to_radians();
            }
            assert!(transform(&src, &p, &mut pt).is_err(), "{dst}");
        }
    }
}
