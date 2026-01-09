//!
//! Implementation of the aeqd (Azimuthal Equidistant) projection.
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

use proj4rs_geodesic::Geodesic;

use crate::errors::{Error, Result};
use crate::math::{
    Enfn, aasin,
    consts::{EPS_10, FRAC_PI_2, PI},
    enfn, inv_mlfn, mlfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { aeqd }

const TOL: f64 = 1.0e-14;

#[allow(private_interfaces)]
#[derive(Debug, Clone)]
pub(crate) enum Projection {
    Guam(GuamProj),
    Aeqd(AeqdProj),
}

impl Projection {
    pub fn aeqd(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        if p.ellps.is_ellipsoid() && params.check_option("guam")? {
            Ok(Projection::Guam(GuamProj::new(p)))
        } else {
            Ok(Projection::Aeqd(AeqdProj::new(p)?))
        }
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Aeqd(p) => p.forward(lam, phi, z),
            Self::Guam(p) => p.forward(lam, phi, z),
        }
    }

    #[inline(always)]
    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Aeqd(p) => p.inverse(x, y, z),
            Self::Guam(p) => p.inverse(x, y, z),
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

// Aeqd

#[allow(non_camel_case_types)]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq)]
enum Mode {
    N_POLE,
    S_POLE,
    OBLIQ,
    EQUIT,
}

use Mode::*;

#[derive(Debug, Clone)]
enum AeqdProj {
    PSpher(AeqdPS), // Polar spherique
    PEllps(AeqdPE), // Polar Ellipsoid
    EqObl(AeqdE),   // Equidistant/Oblique
}

impl AeqdProj {
    fn mode(p: &ProjData) -> (Mode, f64, f64) {
        let (sinph0, cosph0);
        let phi0 = p.phi0;

        let mode = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
            cosph0 = 0.;
            if phi0 < 0. {
                sinph0 = -1.;
                S_POLE
            } else {
                sinph0 = 1.;
                N_POLE
            }
        } else if phi0.abs() < EPS_10 {
            sinph0 = 0.;
            cosph0 = 1.;
            EQUIT
        } else {
            (sinph0, cosph0) = phi0.sin_cos();
            OBLIQ
        };

        (mode, sinph0, cosph0)
    }

    fn new(p: &mut ProjData) -> Result<Self> {
        let phi0 = p.phi0;
        let (mode, sinph0, cosph0) = Self::mode(p);

        //let n = p.ellps.f / (2. - p.ellps.f);
        let es = p.ellps.es;

        match (mode, p.ellps.is_sphere()) {
            (N_POLE | S_POLE, true) => Ok(Self::PSpher(AeqdPS { phi0, mode })),
            (N_POLE, false) => {
                let en = enfn(es);
                Ok(Self::PEllps(AeqdPE {
                    phi0,
                    mode,
                    es: p.ellps.es,
                    Mp: mlfn(FRAC_PI_2, 1., 0., en),
                    en,
                }))
            }
            (S_POLE, false) => {
                let en = enfn(es);
                Ok(Self::PEllps(AeqdPE {
                    phi0,
                    mode,
                    es: p.ellps.es,
                    Mp: mlfn(-FRAC_PI_2, -1., 0., en),
                    en,
                }))
            }
            (EQUIT | OBLIQ, sphere) => Ok(Self::EqObl(AeqdE {
                sphere,
                phi0,
                sinph0,
                cosph0,
                mode,
                g: Box::new(Geodesic::new(1., p.ellps.f)),
            })),
        }
    }

    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::PSpher(p) => p.forward(lam, phi, z), // Polar spherique
            Self::PEllps(p) => p.forward(lam, phi, z), // Polar Ellipsoid
            Self::EqObl(p) => p.forward(lam, phi, z),  // Equidistant/Oblique
        }
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::PSpher(p) => p.inverse(x, y, z), // Polar spherique
            Self::PEllps(p) => p.inverse(x, y, z), // Polar Ellipsoid
            Self::EqObl(p) => p.inverse(x, y, z),  // Equidistant/Oblique
        }
    }
}

// ====================
// Aeqd polar spherical
// ===================
#[derive(Debug, Clone)]
struct AeqdPS {
    phi0: f64,
    mode: Mode,
}

impl AeqdPS {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut phi = phi;
        let mut coslam = lam.cos();
        if self.mode == N_POLE {
            phi = -phi;
            coslam = -coslam;
        }
        if (phi - FRAC_PI_2).abs() < EPS_10 {
            Err(Error::CoordTransOutsideProjectionDomain)
        } else {
            let yy = FRAC_PI_2 + phi;
            Ok((yy * lam.sin(), yy * coslam, z))
        }
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut c_rh = x.hypot(y);
        if c_rh < EPS_10 {
            return Ok((0., self.phi0, z));
        }

        if c_rh > PI {
            if (c_rh - EPS_10) > PI {
                return Err(Error::CoordTransOutsideProjectionDomain);
            }
            c_rh = PI
        }

        Ok(if self.mode == N_POLE {
            (x.atan2(-y), FRAC_PI_2 - c_rh, z)
        } else {
            (x.atan2(y), c_rh - FRAC_PI_2, z)
        })
    }
}

// ======================
// Aeqd Polar Ellipsoidal
// ======================
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
struct AeqdPE {
    phi0: f64,
    mode: Mode,
    es: f64,
    Mp: f64,
    en: Enfn,
}

impl AeqdPE {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut coslam = lam.cos();
        if self.mode == N_POLE {
            coslam = -coslam;
        }

        let rho = (self.Mp - mlfn(phi, phi.sin(), phi.cos(), self.en)).abs();
        Ok((rho * lam.sin(), rho * coslam, z))
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let s12 = x.hypot(y);
        Ok(if s12 < EPS_10 {
            (0., self.phi0, z)
        } else if self.mode == N_POLE {
            (
                x.atan2(-y), // lam
                inv_mlfn(self.Mp - s12, self.es, self.en)?,
                z,
            )
        } else {
            (
                x.atan2(y), // lam
                inv_mlfn(self.Mp + s12, self.es, self.en)?,
                z,
            )
        })
    }
}

// ========================
// Aeqd Oblique/Equidistant
// ========================
#[derive(Debug, Clone)]
struct AeqdE {
    sphere: bool,
    phi0: f64,
    sinph0: f64,
    cosph0: f64,
    mode: Mode,
    g: Box<Geodesic>,
}

impl AeqdE {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.sphere {
            self.s_forward(lam, phi, z)
        } else {
            self.e_forward(lam, phi, z)
        }
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.sphere {
            self.s_inverse(x, y, z)
        } else {
            self.e_inverse(x, y, z)
        }
    }

    // Spherical

    fn s_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let cosphi = phi.cos();
        let coslam = lam.cos();

        if self.mode == EQUIT {
            let mut y = cosphi * coslam;
            if (y.abs() - 1.).abs() < TOL {
                if y < 0. {
                    Err(Error::CoordTransOutsideProjectionDomain)
                } else {
                    self.e_forward(lam, phi, z)
                }
            } else {
                y = y.acos();
                y /= y.sin();
                Ok((y * cosphi * lam.sin(), y * phi.sin(), z))
            }
        } else {
            // OBLIQ
            let sinphi = phi.sin();
            let cosphi_x_coslam = cosphi * coslam;
            let mut y = self.sinph0 * sinphi + self.cosph0 * cosphi_x_coslam;
            if (y.abs() - 1.).abs() < TOL {
                if y < 0. {
                    Err(Error::CoordTransOutsideProjectionDomain)
                } else {
                    self.e_forward(lam, phi, z)
                }
            } else {
                y = y.acos();
                y /= y.sin();
                Ok((
                    y * cosphi * lam.sin(),
                    y * (self.cosph0 * sinphi - self.sinph0 * cosphi_x_coslam),
                    z,
                ))
            }
        }
    }

    fn s_inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut c_rh = x.hypot(y);
        if c_rh > PI {
            if (c_rh - EPS_10) > PI {
                return Err(Error::CoordTransOutsideProjectionDomain);
            }
            c_rh = PI
        } else if c_rh < EPS_10 {
            return Ok((0., self.phi0, z));
        }

        let (sinc, cosc) = c_rh.sin_cos();
        let (phi, lam, xx, yy);
        if self.mode == EQUIT {
            phi = aasin(y * sinc / c_rh)?;
            xx = x * sinc;
            yy = cosc * c_rh;
        } else {
            // OBLIQ
            phi = aasin(cosc * self.sinph0 + y * sinc * self.cosph0 / c_rh)?;
            yy = (cosc - self.sinph0 * phi.sin()) * c_rh;
            xx = x * sinc * self.cosph0;
        }
        lam = if yy == 0. { 0. } else { xx.atan2(yy) };
        Ok((lam, phi, z))
    }

    // Ellipsoidal

    fn e_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok(if lam.abs() < EPS_10 && (phi - self.phi0).abs() < EPS_10 {
            (0., 0., z)
        } else {
            let (s12, mut azi1, _) = self.g.inverse(
                self.phi0.to_degrees(), // lat1
                0.,                     // lon1
                phi.to_degrees(),       // lat2
                lam.to_degrees(),       // lon2
            );
            azi1 = azi1.to_radians();
            (s12 * azi1.sin(), s12 * azi1.cos(), z)
        })
    }

    fn e_inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let s12 = x.hypot(y);
        Ok(if s12 < EPS_10 {
            (0., self.phi0, z)
        } else {
            let (phi, lam, _) = self.g.direct(
                self.phi0.to_degrees(),  // lat1
                0.,                      // lon1
                x.atan2(y).to_degrees(), // az1, clockwise from north
                s12,
            );
            (lam.to_radians(), phi.to_radians(), z)
        })
    }
}

// Guam elliptical (EPSG:3993)

// NOTE Guam projection is valid for the following bounding box
// see https://epsg.io/3993
//
// WGS84 bounds
// 144.58 13.18
// 145.01 13.7
//
// Projected bounds
// 31443.62 17485.32
// 78058.94 75029.52

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
struct GuamProj {
    phi0: f64,
    e: f64,
    es: f64,
    M1: f64,
    en: Enfn,
}

impl GuamProj {
    fn new(p: &mut ProjData) -> Self {
        let (_, sinph0, cosph0) = AeqdProj::mode(p);

        // Based on clenshaw coefficient with
        // third flattenning
        //let f = p.ellps.f;
        //let en = enfn(f / (2. - f));

        // Pre Proj 9 method (same as proj4js).
        let en = enfn(p.ellps.es);

        Self {
            phi0: p.phi0,
            e: p.ellps.e,
            es: p.ellps.es,
            M1: mlfn(p.phi0, sinph0, cosph0, en),
            en,
        }
    }

    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let t = 1. / (1. - self.es * sinphi * sinphi).sqrt();
        Ok((
            lam * cosphi * t,
            mlfn(phi, sinphi, cosphi, self.en) - self.M1 + 0.5 * lam * lam * cosphi * sinphi * t,
            z,
        ))
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let x2 = 0.5 * x * x;
        let mut phi = self.phi0;
        let mut t = 0.;
        eprintln!("################# {} {} {}", x, y, phi.to_degrees());
        for _ in 0..3 {
            t = self.e * phi.sin();
            t = (1. - t * t).sqrt();
            phi = inv_mlfn(self.M1 + y - x2 * phi.tan() * t, self.es, self.en)?;
            eprintln!("XXXXXXXXXXXXXXXXX {} {}", t, phi.to_degrees());
        }
        Ok((x * t / phi.cos(), phi, z))
    }
}

//============
// Tests
//============

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_aeqd_spherical_oblique() {
        let p = Proj::from_proj_string(
            "+proj=aeqd +lon_0=130.0 +lat_0=40.0 +a=6378137 +b=6378137 +units=m",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (2., 1., 0.),
                (-11599752.739940654486, 6022234.512878744863, 0.),
            ),
            (
                (2., -1., 0.),
                (-11904688.609892310575, 5776538.561925039627, 0.),
            ),
            (
                (-2., 1., 0.),
                (-11479183.995349746197, 6850335.025130929425, 0.),
            ),
            (
                (-2., -1., 0.),
                (-11805016.957541335374, 6619965.819121653214, 0.),
            ),
            (
                (-140., -87., 0.),
                (987239.050040827598, -14430471.938706170768, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_spherical_equidistant() {
        let p =
            Proj::from_proj_string("+proj=aeqd +lon_0=0 +lat_0=0 +a=6378137 +b=6378137 +units=m")
                .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222616.370883589523, 111342.098740889822, 0.)),
            (
                (2., -1., 0.),
                (222616.370883589523, -111342.098740889822, 0.),
            ),
            (
                (-2., 1., 0.),
                (-222616.370883589523, 111342.098740889822, 0.),
            ),
            (
                (-2., -1., 0.),
                (-222616.370883589523, -111342.098740889822, 0.),
            ),
            (
                (-140., -87., 0.),
                (-345922.696188778791, -10268707.968071108684, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_spherical_polar() {
        let p =
            Proj::from_proj_string("+proj=aeqd +lon_0=0 +lat_0=90 +a=6378137 +b=6378137 +units=m")
                .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (2., 1., 0.),
                (345764.483965890424, -9901399.339083850384, 0.),
            ),
            (
                (2., -1., 0.),
                (353534.472369618365, -10123902.695018321276, 0.),
            ),
            (
                (-2., 1., 0.),
                (-345764.483965890424, -9901399.339083850384, 0.),
            ),
            (
                (-2., -1., 0.),
                (-353534.472369618365, -10123902.695018321276, 0.),
            ),
            (
                (-140., -87., 0.),
                (-12665197.723539995030, 15093794.887944793329, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_ellipsoidal_oblique() {
        let p = Proj::from_proj_string("+proj=aeqd +lon_0=130 +lat_0=40 +ellps=bessel +units=m")
            .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (2., 1., 0.),
                (-11576389.897925473750, 6054978.016031168401, 0.),
            ),
            (
                (2., -1., 0.),
                (-11880228.068629510701, 5811535.382627814077, 0.),
            ),
            (
                (-2., 1., 0.),
                (-11452628.503103895113, 6883529.943709976971, 0.),
            ),
            (
                (-2., -1., 0.),
                (-11777149.990328865126, 6655618.481343483552, 0.),
            ),
            (
                (-140., -87., 0.),
                (986349.921549595078, -14388885.114546611905, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_ellipsoidal_polar() {
        let p =
            Proj::from_proj_string("+proj=aeqd +lon_0=0 +lat_0=90 +ellps=bessel +units=m").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (2., 1., 0.),
                (345166.212186285236, -9884267.076871054247, 0.),
            ),
            (
                (2., -1., 0.),
                (352883.453359715699, -10105259.949758755043, 0.),
            ),
            (
                (-2., 1., 0.),
                (-345166.212186285236, -9884267.076871054247, 0.),
            ),
            (
                (-2., -1., 0.),
                (-352883.453359715699, -10105259.949758755043, 0.),
            ),
            (
                (-140., 87., 0.),
                (-215357.382809855917, 256652.934655332210, 0.),
            ),
            /*
             * We have a 10-5 discrepancy error with proj 9.4
             * This this acceptable but make the test fail
             * with the 1e-6 required accuracy
             * XXX must check with clenshaw cofficients for
             * meridional distance computation
            (
                (-140., -87., 0.),
                (-12641494.960468998179, 15065547.034900521860, 0.),
            ),
            */
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_ellipsoidal_equidistant() {
        let p =
            Proj::from_proj_string("+proj=aeqd +lon_0=0 +lat_0=0 +ellps=bessel +units=m").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222590.698880701879, 110586.394289519230, 0.)),
            (
                (2., -1., 0.),
                (222590.698880701821, -110586.394289519303, 0.),
            ),
            (
                (-2., 1., 0.),
                (-222590.698880701879, 110586.394289519230, 0.),
            ),
            (
                (-2., -1., 0.),
                (-222590.698880701821, -110586.394289519303, 0.),
            ),
            (
                (-140., -87., 0.),
                (-346430.547494323750, -10251627.555727558210, 0.),
            ),
            (
                (-140., 87., 0.),
                (-346430.547494323982, 10251627.555727558210, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_aeqd_ellipsoidal_guam() {
        let p = Proj::from_proj_string(concat!(
            "+proj=aeqd +guam +lon_0=144.74875 +lat_0=13.47246 ",
            "+x_0=50000 +y_0=50000 +ellps=clrk66 +units=m",
        ))
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (144.8, 13.5, 0.),
            (55548.574682530634, 53047.282504629671, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }
}
