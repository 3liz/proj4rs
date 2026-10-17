//!
//! Orthographic
//!
//! ref: <https://proj.org/en/stable/operations/projections/ortho.html>
//!
//! ortho: "Orthographic" "\n\tAzi, Sph&Ell";
//!
//! Ported from PROJ `src/projections/ortho.cpp`. The ellipsoidal form follows
//! EPSG guidance note 7-2 (March 2020), section 3.3.5, as PROJ does since 7.2.
//!
//! Parameters:
//!
//! * lat_0, lon_0, x_0, y_0, k_0
//! * alpha: azimuth of the projection plane's y axis (rotation, degrees)
//!
//! Copyright (c) 2000, Frank Warmerdam
//! Copyright (c) 2020, Even Rouault
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
    adjlon,
    consts::{EPS_10, FRAC_PI_2, PI},
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { ortho }

const MAX_ITER: usize = 20;

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
    k0: f64,
    sinalpha: f64,
    cosalpha: f64,
    mode: Mode,
    // Ellipsoidal parameters; `es == 0` selects the spherical form
    es: f64,
    one_es: f64,
    b_a: f64,
    nu0: f64,
    y_shift: f64,
    y_scale: f64,
}

impl Projection {
    pub fn ortho(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let phi0 = p.phi0;
        let (sinph0, cosph0) = phi0.sin_cos();

        let mode = if (phi0.abs() - FRAC_PI_2).abs() <= EPS_10 {
            if phi0 < 0. {
                Mode::S_POLE
            } else {
                Mode::N_POLE
            }
        } else if phi0.abs() > EPS_10 {
            Mode::OBLIQ
        } else {
            Mode::EQUIT
        };

        let (es, one_es) = (p.ellps.es, p.ellps.one_es);
        let (nu0, y_shift, y_scale) = if es == 0. {
            (1., 0., 1.)
        } else {
            let nu0 = 1. / (1. - es * sinph0 * sinph0).sqrt();
            (
                nu0,
                es * nu0 * sinph0 * cosph0,
                1. / (1. - es * cosph0 * cosph0).sqrt(),
            )
        };

        let alpha = params.try_angular_value("alpha")?.unwrap_or(0.);
        let (sinalpha, cosalpha) = alpha.sin_cos();

        Ok(Self {
            phi0,
            sinph0,
            cosph0,
            k0: p.k0,
            sinalpha,
            cosalpha,
            mode,
            es,
            one_es,
            b_a: p.ellps.b / p.ellps.a,
            nu0,
            y_shift,
            y_scale,
        })
    }

    #[inline(always)]
    fn rotate(&self, xp: f64, yp: f64) -> (f64, f64) {
        (
            (xp * self.cosalpha - yp * self.sinalpha) * self.k0,
            (xp * self.sinalpha + yp * self.cosalpha) * self.k0,
        )
    }

    #[inline(always)]
    fn unrotate(&self, x: f64, y: f64) -> (f64, f64) {
        (
            (self.cosalpha * x + self.sinalpha * y) / self.k0,
            (-self.sinalpha * x + self.cosalpha * y) / self.k0,
        )
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.es == 0. {
            self.s_forward(lam, phi, z)
        } else {
            self.e_forward(lam, phi, z)
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (x, y) = self.unrotate(x, y);
        let (lam, phi) = if self.es == 0. {
            self.s_inverse(x, y)?
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

    fn s_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let mut coslam = lam.cos();

        let yp = match self.mode {
            Mode::EQUIT => {
                if cosphi * coslam < -EPS_10 {
                    return Err(Error::CoordinateOutOfRange);
                }
                sinphi
            }
            Mode::OBLIQ => {
                // Visibility: dot product of the normals at the centre and
                // at the point (Snyder, eq. 5-3).
                if self.sinph0 * sinphi + self.cosph0 * cosphi * coslam < -EPS_10 {
                    return Err(Error::CoordinateOutOfRange);
                }
                self.cosph0 * sinphi - self.sinph0 * cosphi * coslam
            }
            Mode::N_POLE | Mode::S_POLE => {
                if self.mode == Mode::N_POLE {
                    coslam = -coslam;
                }
                if (phi - self.phi0).abs() - EPS_10 > FRAC_PI_2 {
                    return Err(Error::CoordinateOutOfRange);
                }
                cosphi * coslam
            }
        };
        let xp = cosphi * lam.sin();
        let (x, y) = self.rotate(xp, yp);
        Ok((x, y, z))
    }

    fn s_inverse(&self, mut x: f64, mut y: f64) -> Result<(f64, f64)> {
        let rh = x.hypot(y);
        let mut sinc = rh;
        if sinc > 1. {
            if (sinc - 1.) > EPS_10 {
                return Err(Error::CoordinateOutOfRange);
            }
            sinc = 1.;
        }
        let cosc = (1. - sinc * sinc).sqrt();

        if rh.abs() <= EPS_10 {
            return Ok((0., self.phi0));
        }

        let phi = match self.mode {
            Mode::N_POLE => {
                y = -y;
                sinc.acos()
            }
            Mode::S_POLE => -sinc.acos(),
            Mode::EQUIT | Mode::OBLIQ => {
                let phi = if self.mode == Mode::EQUIT {
                    let phi = y * sinc / rh;
                    x *= sinc;
                    y = cosc * rh;
                    phi
                } else {
                    let phi = cosc * self.sinph0 + y * sinc * self.cosph0 / rh;
                    y = (cosc - self.sinph0 * phi) * rh;
                    x *= sinc * self.cosph0;
                    phi
                };
                if phi.abs() >= 1. {
                    FRAC_PI_2.copysign(phi)
                } else {
                    phi.asin()
                }
            }
        };

        let lam = if y == 0. && matches!(self.mode, Mode::OBLIQ | Mode::EQUIT) {
            if x == 0. {
                0.
            } else if x < 0. {
                -FRAC_PI_2
            } else {
                FRAC_PI_2
            }
        } else {
            x.atan2(y)
        };
        Ok((lam, phi))
    }

    //
    // Ellipsoid
    //

    #[inline(always)]
    fn e_xy(&self, sinphi: f64, cosphi: f64, sinlam: f64, coslam: f64) -> (f64, f64) {
        let nu = 1. / (1. - self.es * sinphi * sinphi).sqrt();
        (
            nu * cosphi * sinlam,
            nu * (sinphi * self.cosph0 - cosphi * self.sinph0 * coslam)
                + self.es * (self.nu0 * self.sinph0 - nu * sinphi) * self.cosph0,
        )
    }

    fn e_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let (sinlam, coslam) = lam.sin_cos();

        // Same visibility condition as the spherical case
        if self.sinph0 * sinphi + self.cosph0 * cosphi * coslam < -EPS_10 {
            return Err(Error::CoordinateOutOfRange);
        }
        let (xp, yp) = self.e_xy(sinphi, cosphi, sinlam, coslam);
        let (x, y) = self.rotate(xp, yp);
        Ok((x, y, z))
    }

    fn e_inverse(&self, x: f64, y: f64) -> Result<(f64, f64)> {
        let es = self.es;

        match self.mode {
            Mode::N_POLE | Mode::S_POLE => {
                // rh^2 = cosphi^2 / (1 - es sinphi^2)
                // => cosphi^2 = rh^2 (1 - es) / (1 - es rh^2)
                let sign = if self.mode == Mode::N_POLE { 1. } else { -1. };
                let rh2 = x * x + y * y;
                let phi = if rh2 >= 1. - 1e-15 {
                    if (rh2 - 1.) > EPS_10 {
                        return Err(Error::CoordinateOutOfRange);
                    }
                    0.
                } else {
                    (rh2 * self.one_es / (1. - es * rh2)).sqrt().acos() * sign
                };
                Ok((x.atan2(-y * sign), phi))
            }
            Mode::EQUIT => {
                // Equation of the ellipse
                let yr = y / self.b_a;
                if x * x + yr * yr > 1. + 1e-11 {
                    return Err(Error::CoordinateOutOfRange);
                }
                let sinphi2 = if y == 0. {
                    0.
                } else {
                    let t = self.one_es / y;
                    1. / (t * t + es)
                };
                if sinphi2 > 1. - 1e-11 {
                    return Ok((0., FRAC_PI_2.copysign(y)));
                }
                let phi = sinphi2.sqrt().asin().copysign(y);
                let sinlam = x * ((1. - es * sinphi2) / (1. - sinphi2)).sqrt();
                let lam = if sinlam.abs() - 1. > -1e-15 {
                    FRAC_PI_2.copysign(x)
                } else {
                    sinlam.asin()
                };
                Ok((lam, phi))
            }
            Mode::OBLIQ => {
                let yr = (y - self.y_shift) / self.y_scale;
                if x * x + yr * yr > 1. + 1e-11 {
                    return Err(Error::CoordinateOutOfRange);
                }

                // Newton iteration (EPSG guidance note 7-2, 3.3.5), seeded
                // with the spherical inverse, which converges near the poles.
                let (mut lam, mut phi) = self.s_inverse(x, yr)?;

                for _ in 0..MAX_ITER {
                    let (sinphi, cosphi) = phi.sin_cos();
                    let (sinlam, coslam) = lam.sin_cos();
                    let one_minus_es_sinphi2 = 1. - es * sinphi * sinphi;
                    let nu = 1. / one_minus_es_sinphi2.sqrt();
                    let (xn, yn) = self.e_xy(sinphi, cosphi, sinlam, coslam);
                    let rho = self.one_es * nu / one_minus_es_sinphi2;
                    let j11 = -rho * sinphi * sinlam;
                    let j12 = nu * cosphi * coslam;
                    let j21 = rho * (cosphi * self.cosph0 + sinphi * self.sinph0 * coslam);
                    let j22 = nu * self.sinph0 * cosphi * sinlam;
                    let d = j11 * j22 - j12 * j21;
                    let (dx, dy) = (x - xn, y - yn);
                    let dphi = (j22 * dx - j12 * dy) / d;
                    let dlam = (-j21 * dx + j11 * dy) / d;
                    phi += dphi;
                    if phi > FRAC_PI_2 {
                        phi = FRAC_PI_2 - (phi - FRAC_PI_2);
                        lam = adjlon(lam + PI);
                    } else if phi < -FRAC_PI_2 {
                        phi = -FRAC_PI_2 + (-FRAC_PI_2 - phi);
                        lam = adjlon(lam + PI);
                    }
                    lam += dlam;
                    if dphi.abs() < 1e-12 && dlam.abs() < 1e-12 {
                        return Ok((lam, phi));
                    }
                }
                Err(Error::InverseProjectionFailure)
            }
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
    fn proj_ortho_s_obliq() {
        let p = Proj::from_proj_string("+proj=ortho +R=6370997 +lat_0=-42 +lon_0=147").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (288760.773147241, 1319545.149494376, 0.)),
            (
                (100., -60., 0.),
                (-2329726.119477315, -2646571.464291994, 0.),
            ),
            (
                (-133., -37., 0.),
                (5010804.699273809, -2258133.734182011, 0.),
            ),
            ((147., 0., 0.), (0., 4263029.085720466, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_s_equit() {
        let p = Proj::from_proj_string("+proj=ortho +R=6370997 +lat_0=0 +lon_0=0").unwrap();

        let inputs = [
            ((10., 20., 0.), (1039593.240514260, 2179009.307067405, 0.)),
            ((-80., 5., 0.), (-6250331.987234295, 555268.975578102, 0.)),
            ((45., -45., 0.), (3185498.5, -4504975.181619151, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_s_poles() {
        let p = Proj::from_proj_string("+proj=ortho +R=6370997 +lat_0=90").unwrap();
        let inputs = [
            ((10., 60., 0.), (553156.009485740, -3137103.620008760, 0.)),
            (
                (-120., 10., 0.),
                (-5433622.858463421, 3137103.620008758, 0.),
            ),
            ((170., 85., 0.), (96421.445724121, 546833.192156461, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string("+proj=ortho +R=6370997 +lat_0=-90 +lon_0=30").unwrap();
        let inputs = [
            ((10., -60., 0.), (-1089504.653533703, 2993389.433974581, 0.)),
            (
                (-120., -10., 0.),
                (-3137103.620008761, -5433622.858463419, 0.),
            ),
            ((170., -85., 0.), (356920.017544941, -425360.713177972, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_e_obliq() {
        let p = Proj::from_proj_string("+proj=ortho +ellps=WGS84 +lat_0=-42 +lon_0=147").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (289326.597573795, 1316745.581434381, 0.)),
            (
                (100., -60., 0.),
                (-2338214.274110457, -2649929.082205364, 0.),
            ),
            (
                (-133., -37., 0.),
                (5022512.762328568, -2265554.666728071, 0.),
            ),
            ((147., 0., 0.), (0., 4246542.902040465, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_e_equit() {
        let p = Proj::from_proj_string("+proj=ortho +ellps=WGS84 +lat_0=0 +lon_0=0").unwrap();

        let inputs = [
            ((10., 20., 0.), (1041166.061652678, 2167696.787828757, 0.)),
            ((-80., 5., 0.), (-6257495.860544777, 552183.960027770, 0.)),
            ((45., -45., 0.), (3194419.145060574, -4487348.408865920, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_e_poles() {
        let p = Proj::from_proj_string("+proj=ortho +ellps=WGS84 +lat_0=90").unwrap();
        let inputs = [
            ((10., 60., 0.), (555171.385329928, -3148533.384393596, 0.)),
            (
                (-120., 10., 0.),
                (-5440261.453779826, 3140936.414801726, 0.),
            ),
            ((170., 85., 0.), (96851.760438474, 549273.628172689, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string("+proj=ortho +ellps=WGS84 +lat_0=-90 +lon_0=30").unwrap();
        let inputs = [
            ((10., -60., 0.), (-1093474.169046881, 3004295.588213213, 0.)),
            (
                (-120., -10., 0.),
                (-3140936.414801729, -5440261.453779825, 0.),
            ),
            ((170., -85., 0.), (358512.899027308, -427259.035407160, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_e_alpha_k0() {
        let p = Proj::from_proj_string(
            "+proj=ortho +ellps=WGS84 +lat_0=40 +lon_0=-100 +alpha=30 +k_0=0.9996 +x_0=1000 +y_0=-500",
        )
        .unwrap();

        let inputs = [
            ((-100., 40., 0.), (1000., -500., 0.)),
            ((-90., 45., 0.), (380806.930194124, 909970.369811893, 0.)),
            (
                (-120., 20., 0.),
                (-804625.158396952, -2704869.375100226, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_ortho_hidden_hemisphere() {
        use crate::transform::transform;

        let src = Proj::from_proj_string("+proj=longlat +ellps=WGS84").unwrap();
        for def in [
            "+proj=ortho +R=6370997 +lat_0=-42 +lon_0=147",
            "+proj=ortho +ellps=WGS84 +lat_0=-42 +lon_0=147",
            "+proj=ortho +ellps=WGS84 +lat_0=90",
        ] {
            let dst = Proj::from_proj_string(def).unwrap();
            // Antipode of the centre, or a point on the far hemisphere
            let mut pt = ((-33f64).to_radians(), (42f64).to_radians(), 0.);
            if def.contains("lat_0=90") {
                pt.1 = (-10f64).to_radians();
            }
            assert!(transform(&src, &dst, &mut pt).is_err(), "{def}");
        }
    }
}
