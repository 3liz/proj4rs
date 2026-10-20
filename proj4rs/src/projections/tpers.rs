//!
//! Tilted perspective and near-sided perspective
//!
//! ref: <https://proj.org/en/stable/operations/projections/tpers.html>
//! ref: <https://proj.org/en/stable/operations/projections/nsper.html>
//!
//! tpers: "Tilted perspective" "\n\tAzi, Sph\n\ttilt= azi= h=";
//! nsper: "Near-sided perspective" "\n\tAzi, Sph\n\th=";
//!
//! Ported from PROJ `src/projections/nsper.cpp`. Spherical only: on an
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
use crate::math::consts::{EPS_10, FRAC_PI_2};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { tpers, nsper }

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
    p: f64,
    rp: f64,
    pn1: f64,
    pfact: f64,
    h: f64,
    cg: f64,
    sg: f64,
    sw: f64,
    cw: f64,
    mode: Mode,
    tilt: bool,
}

impl Projection {
    pub fn tpers(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let omega = params.try_angular_value("tilt")?.unwrap_or(0.);
        let gamma = params.try_angular_value("azi")?.unwrap_or(0.);
        Self::init(p, params, true, omega, gamma)
    }

    pub fn nsper(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::init(p, params, false, 0., 0.)
    }

    fn init(p: &ProjData, params: &ParamList, tilt: bool, omega: f64, gamma: f64) -> Result<Self> {
        let height: f64 = params.try_value("h")?.unwrap_or(0.);

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

        // Normalize by radius
        let pn1 = height / p.ellps.a;
        if pn1 <= 0. || pn1 > 1e10 {
            return Err(Error::InvalidParameterValue("Invalid value for 'h'."));
        }
        let pp = 1. + pn1;
        let h = 1. / pn1;
        let (sg, cg) = gamma.sin_cos();
        let (sw, cw) = omega.sin_cos();

        Ok(Self {
            phi0,
            sinph0,
            cosph0,
            p: pp,
            rp: 1. / pp,
            pn1,
            pfact: (pp + 1.) * h,
            h,
            cg,
            sg,
            sw,
            cw,
            mode,
            tilt,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let coslam = lam.cos();

        let cosz = match self.mode {
            Mode::OBLIQ => self.sinph0 * sinphi + self.cosph0 * cosphi * coslam,
            Mode::EQUIT => cosphi * coslam,
            Mode::S_POLE => -sinphi,
            Mode::N_POLE => sinphi,
        };
        if cosz < self.rp {
            return Err(Error::CoordinateOutOfRange);
        }

        let r = self.pn1 / (self.p - cosz);
        let mut x = r * cosphi * lam.sin();
        let mut y = r * match self.mode {
            Mode::OBLIQ => self.cosph0 * sinphi - self.sinph0 * cosphi * coslam,
            Mode::EQUIT => sinphi,
            Mode::N_POLE => -cosphi * coslam,
            Mode::S_POLE => cosphi * coslam,
        };

        if self.tilt {
            let yt = y * self.cg + x * self.sg;
            let ba = 1. / (yt * self.sw * self.h + self.cw);
            x = (x * self.cg - y * self.sg) * self.cw * ba;
            y = yt * ba;
        }
        Ok((x, y, z))
    }

    pub fn inverse(&self, mut x: f64, mut y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.tilt {
            let yt = 1. / (self.pn1 - y * self.sw);
            let bm = self.pn1 * x * yt;
            let bq = self.pn1 * y * self.cw * yt;
            x = bm * self.cg + bq * self.sg;
            y = bq * self.cg - bm * self.sg;
        }

        let rh = x.hypot(y);
        if rh.abs() <= EPS_10 {
            return Ok((0., self.phi0, z));
        }

        let sinz = 1. - rh * rh * self.pfact;
        if sinz < 0. {
            return Err(Error::CoordinateOutOfRange);
        }
        let sinz = (self.p - sinz.sqrt()) / (self.pn1 / rh + rh / self.pn1);
        let cosz = (1. - sinz * sinz).sqrt();

        let phi = match self.mode {
            Mode::OBLIQ => {
                let phi = (cosz * self.sinph0 + y * sinz * self.cosph0 / rh).asin();
                y = (cosz - self.sinph0 * phi.sin()) * rh;
                x *= sinz * self.cosph0;
                phi
            }
            Mode::EQUIT => {
                let phi = (y * sinz / rh).asin();
                y = cosz * rh;
                x *= sinz;
                phi
            }
            Mode::N_POLE => {
                y = -y;
                cosz.asin()
            }
            Mode::S_POLE => -cosz.asin(),
        };
        Ok((x.atan2(y), phi, z))
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

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_tpers_obliq() {
        let p = Proj::from_proj_string(
            "+proj=tpers +R=6370997 +lat_0=-42 +lon_0=147 +h=3000000 +tilt=30 +azi=45",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (-575211.027126064, 1036328.413044423, 0.)),
            ((140., -50., 0.), (344750.161414741, -1371938.406334842, 0.)),
            ((160., -40., 0.), (564968.408922931, 840956.192087849, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_tpers_poles() {
        let p = Proj::from_proj_string(
            "+proj=tpers +ellps=WGS84 +lat_0=90 +h=1000000 +tilt=20 +azi=-30",
        )
        .unwrap();
        let inputs = [
            ((10., 80., 0.), (-527515.170085313, -1542351.175448220, 0.)),
            ((-120., 75., 0.), (-474999.365277009, 875523.565923786, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string(
            "+proj=tpers +R=6370997 +lat_0=-90 +lon_0=30 +h=5000000 +tilt=10 +azi=100",
        )
        .unwrap();
        let inputs = [
            (
                (10., -60., 0.),
                (-2475208.711844577, -1451108.007337406, 0.),
            ),
            (
                (-120., -40., 0.),
                (3284507.482964687, -1213904.901100970, 0.),
            ),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_tpers_equit() {
        let p = Proj::from_proj_string("+proj=tpers +ellps=WGS84 +lat_0=0 +lon_0=0 +h=35786000")
            .unwrap();
        let inputs = [
            ((10., 20., 0.), (1027105.010111250, 2152833.713367413, 0.)),
            (
                (-60., -30., 0.),
                (-4344566.422358030, -2896377.614905353, 0.),
            ),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_nsper() {
        let p = Proj::from_proj_string("+proj=nsper +ellps=WGS84 +lat_0=40 +lon_0=-100 +h=3000000")
            .unwrap();
        let inputs = [
            ((-100., 40., 0.), (0., 0., 0.)),
            ((-90., 45., 0.), (763619.441029603, 584966.429458347, 0.)),
            (
                (-120., 25., 0.),
                (-1702224.772677926, -1228368.435913821, 0.),
            ),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string("+proj=nsper +R=6370997 +lat_0=90 +h=2000000").unwrap();
        let inputs = [
            ((10., 60., 0.), (387696.497444708, -2198736.097487637, 0.)),
            (
                (-120., 50., 0.),
                (-2032095.578616672, 1173230.929333383, 0.),
            ),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_tpers_hidden_hemisphere() {
        use crate::transform::transform;

        for (src, dst) in [
            (
                "+proj=longlat +R=6370997",
                "+proj=tpers +R=6370997 +lat_0=-42 +lon_0=147 +h=3000000 +tilt=30 +azi=45",
            ),
            (
                "+proj=longlat +ellps=WGS84",
                "+proj=nsper +ellps=WGS84 +lat_0=40 +lon_0=-100 +h=3000000",
            ),
        ] {
            let src = Proj::from_proj_string(src).unwrap();
            let p = Proj::from_proj_string(dst).unwrap();
            let mut pt = ((-33f64).to_radians(), (42f64).to_radians(), 0.);
            assert!(transform(&src, &p, &mut pt).is_err(), "{dst}");
        }
    }

    #[test]
    fn proj_tpers_invalid_h() {
        for def in [
            "+proj=tpers +ellps=WGS84 +lat_0=40",
            "+proj=nsper +ellps=WGS84 +lat_0=40 +h=-10",
        ] {
            assert!(Proj::from_proj_string(def).is_err(), "{def}");
        }
    }
}
