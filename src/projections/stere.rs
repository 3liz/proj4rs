//!
//! From proj/stere.c
//!
//! see also https://proj.org/operations/projections/stere.html
//!
//! Stereographic - Azimuthal
//!
//! Provide forward and inverse spherical and ellipsoidal projection.
//! Defined area: global 2D.
//!
//! Geodetique coordinates to projected coordinates.
//!
//! Universal Polar Stereographic (Azimuthal, Spherical, ellipsoidal)
//!
//! ### Derived projection:
//!
//! ups:  Universal Polar Stereographic
//!
//! see also https://proj.org/operations/projections/ups.html
//!
//! stere: "Stereographic" "\n\tAzi, Sph&Ell\n\tlat_ts=";
//! ups: "Universal Polar Stereographic") "\n\tAzi, Sph&Ell\n\tsouth";
//!
use crate::errors::{Error, Result};
use crate::math::{
    consts::{EPS_10, FRAC_PI_2, FRAC_PI_4},
    tsfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

#[inline]
fn ssfn(phit: f64, sinphi: f64, eccen: f64) -> f64 {
    let sinphi = sinphi * eccen;
    (0.5 * (FRAC_PI_2 + phit)).tan() * ((1. - sinphi) / (1. + sinphi)).powf(0.5 * eccen)
}

#[allow(non_camel_case_types)]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Copy, Clone, PartialEq)]
enum Mode {
    S_POLE,
    N_POLE,
    OBLIQ,
    EQUIT,
}

use Mode::*;

// Projection stub
super::projection! { stere, ups }

#[derive(Debug)]
pub(crate) struct Projection {
    mode: Mode,
    e: f64,
    phi0: f64,
    sinx1: f64,
    cosx1: f64,
    akm1: f64,
}

impl Projection {
    #[inline]
    pub fn is_ellipsoid(&self) -> bool {
        self.e != 0.
    }

    // -----------
    // stere
    // -----------
    pub fn stere(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::init(p, params.try_value("lat_ts")?.unwrap_or(FRAC_PI_2))
    }

    // -----------
    // ups
    // -----------
    pub fn ups(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        // International Ellipsoid
        p.phi0 = if params.check_option("south")? {
            -FRAC_PI_2
        } else {
            FRAC_PI_2
        };

        if p.ellps.is_sphere() {
            Err(Error::EllipsoidRequired)
        } else {
            p.k0 = 0.994;
            p.x0 = 2_000_000.;
            p.y0 = 2_000_000.;
            p.lam0 = 0.;
            Self::init(p, FRAC_PI_2)
        }
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.is_ellipsoid() {
            self.e_forward(lam, phi, z)
        } else {
            self.s_forward(lam, phi, z)
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.is_ellipsoid() {
            self.e_inverse(x, y, z)
        } else {
            self.s_inverse(x, y, z)
        }
    }

    //------------------
    // Ellipsoidal
    //------------------

    #[inline(always)]
    pub fn e_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let coslam = lam.cos();
        let sinlam = lam.sin();
        let sinphi = phi.sin();

        let (x, y) = match self.mode {
            OBLIQ => {
                let x = 2. * ssfn(phi, sinphi, self.e).atan() - FRAC_PI_2;
                let (sinx, cosx) = x.sin_cos();
                let denom = self.cosx1 * (1. + self.sinx1 * sinx + self.cosx1 * cosx * coslam);
                if denom == 0. {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                let a = self.akm1 / denom;
                (
                    a * cosx,
                    a * (self.cosx1 * sinx - self.sinx1 * cosx * coslam),
                )
            }
            EQUIT => {
                let x = 2. * ssfn(phi, sinphi, self.e).atan() - FRAC_PI_2;
                let (sinx, cosx) = x.sin_cos();
                // avoid division by zero
                let denom = 1. + cosx * coslam;
                if denom == 0. {
                    return Err(Error::ToleranceConditionError);
                }
                let a = self.akm1 / denom;
                (a * cosx, a * sinx)
            }
            S_POLE => {
                if (phi.abs() - FRAC_PI_2).abs() < 1e-15 {
                    (0., 0.)
                } else {
                    let x = self.akm1 * tsfn(-phi, -sinphi, self.e);
                    (x, x * coslam)
                }
            }
            N_POLE => {
                if (phi.abs() - FRAC_PI_2).abs() < 1e-15 {
                    (0., 0.)
                } else {
                    let x = self.akm1 * tsfn(phi, sinphi, self.e);
                    (x, -x * coslam)
                }
            }
        };
        Ok((x * sinlam, y, z))
    }

    #[inline(always)]
    pub fn e_inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rho = x.hypot(y);

        let (halfpi, halfe, tp, mut phi_l);

        let (xx, yy) = match self.mode {
            OBLIQ | EQUIT => {
                let (sinphi, cosphi) = (2. * (rho * self.cosx1).atan2(self.akm1)).sin_cos();
                phi_l = if rho == 0.0 {
                    (cosphi * self.sinx1).asin()
                } else {
                    (cosphi * self.sinx1 + (y * sinphi * self.cosx1 / rho)).asin()
                };
                tp = (0.5 * (FRAC_PI_2 + phi_l)).tan();
                halfpi = FRAC_PI_2;
                halfe = 0.5 * self.e;
                (
                    x * sinphi,
                    rho * self.cosx1 * cosphi - y * self.sinx1 * sinphi,
                )
            }
            _ => {
                tp = -rho / self.akm1;
                phi_l = FRAC_PI_2 - 2. * tp.atan();
                halfpi = -FRAC_PI_2;
                halfe = -0.5 * self.e;
                if self.mode == N_POLE {
                    (x, -y)
                } else {
                    (x, y)
                }
            }
        };

        const NITER: usize = 8;

        let (mut lam, mut phi) = (0., 0.);
        let mut i = NITER;
        while i > 0 {
            let sinphi = self.e * phi_l.sin();
            phi = 2. * (tp * ((1. + sinphi) / (1. - sinphi)).powf(halfe)).atan() - halfpi;
            if (phi_l - phi).abs() < EPS_10 {
                if self.mode == S_POLE {
                    phi = -phi
                }
                lam = if xx == 0. && yy == 0. {
                    0.
                } else {
                    xx.atan2(yy)
                };
                break;
            }
            phi_l = phi;
            i -= 1;
        }

        if i == 0 {
            Err(Error::CoordTransOutsideProjectionDomain)
        } else {
            Ok((lam, phi, z))
        }
    }

    //------------------
    // Spherical
    //------------------
    #[inline(always)]
    pub fn s_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let (sinlam, coslam) = lam.sin_cos();

        let mode = self.mode;
        let (x, y) = match mode {
            EQUIT | OBLIQ => {
                let (mut y, fac) = if mode == EQUIT {
                    (1. + cosphi * coslam, sinphi)
                } else {
                    (
                        1. + self.sinx1 * sinphi + self.cosx1 * cosphi * coslam,
                        self.cosx1 * sinphi - self.sinx1 * cosphi * coslam,
                    )
                };
                if y <= EPS_10 {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                y = self.akm1 / y;
                (y * cosphi * sinlam, y * fac)
            }
            _ => {
                let (phi, coslam) = if mode == N_POLE {
                    (-phi, -coslam)
                } else {
                    (phi, coslam)
                };
                if (phi - FRAC_PI_2).abs() < 1.0e-8 {
                    return Err(Error::CoordTransOutsideProjectionDomain);
                }
                let y = self.akm1 * (FRAC_PI_4 + 0.5 * phi).tan();
                (y * sinlam, y * coslam)
            }
        };
        Ok((x, y, z))
    }

    #[inline(always)]
    pub fn s_inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rh = x.hypot(y);

        let (sinc, cosc) = (2. * (rh / self.akm1).atan()).sin_cos();
        let (lam, phi) = match self.mode {
            EQUIT => (
                // lam
                if cosc != 0. || x != 0. {
                    (x * sinc).atan2(cosc * rh)
                } else {
                    0.
                },
                // phi
                if rh.abs() <= EPS_10 {
                    0.
                } else {
                    (y * sinc / rh).asin()
                },
            ),
            OBLIQ => {
                let phi = if rh.abs() <= EPS_10 {
                    self.phi0
                } else {
                    (cosc * self.sinx1 + y * sinc * self.cosx1 / rh).asin()
                };
                let c = cosc - self.sinx1 * phi.sin();
                (
                    // lam
                    if c != 0. || x != 0. {
                        (x * sinc * self.cosx1).atan2(c * rh)
                    } else {
                        0.
                    },
                    phi,
                )
            }
            N_POLE => (
                // lam
                if x == 0. && y == 0. { 0. } else { x.atan2(-y) },
                // phi
                if rh.abs() <= EPS_10 {
                    self.phi0
                } else {
                    cosc.asin()
                },
            ),
            S_POLE => (
                // lam
                if x == 0. && y == 0. { 0. } else { x.atan2(y) },
                // phi
                if rh.abs() <= EPS_10 {
                    self.phi0
                } else {
                    (-cosc).asin()
                },
            ),
        };

        Ok((lam, phi, z))
    }

    //-----------------------
    // General initialization
    //-----------------------
    pub fn init(p: &mut ProjData, phits: f64) -> Result<Self> {
        let t = p.phi0.abs();
        let mode = if (t - FRAC_PI_2).abs() < EPS_10 {
            if p.phi0 < 0. {
                S_POLE
            } else {
                N_POLE
            }
        } else if t > EPS_10 {
            OBLIQ
        } else {
            EQUIT
        };

        let phits = phits.abs();

        let (mut sinx1, mut cosx1) = (0., 0.);

        // Get ellipsoid
        let el = &p.ellps;

        let akm1 = if el.is_ellipsoid() {
            // Ellipsoidal
            let ecc = el.e;
            match mode {
                N_POLE | S_POLE => {
                    if (phits - FRAC_PI_2).abs() < EPS_10 {
                        2. * p.k0 / ((1. + ecc).powf(1. + ecc) * (1. - ecc).powf(1. - ecc)).sqrt()
                    } else {
                        let s = phits.sin();
                        let t = s * ecc;
                        phits.cos() / tsfn(phits, s, ecc) / (1. - t * t).sqrt()
                    }
                }
                EQUIT | OBLIQ => {
                    let mut t = p.phi0.sin();
                    let x = 2. * ssfn(p.phi0, t, ecc).atan() - FRAC_PI_2;
                    sinx1 = x.sin();
                    cosx1 = x.cos();
                    t *= ecc;
                    2. * p.k0 * p.phi0.cos() / (1. - t * t).sqrt()
                }
            }
        } else {
            match mode {
                EQUIT => 2. * p.k0,
                OBLIQ => {
                    sinx1 = p.phi0.sin();
                    cosx1 = p.phi0.cos();
                    2. * p.k0
                }
                S_POLE | N_POLE => {
                    if (phits - FRAC_PI_2).abs() >= EPS_10 {
                        phits.cos() / (FRAC_PI_4 - 0.5 * phits).tan()
                    } else {
                        2. * p.k0
                    }
                }
            }
        };

        Ok(Self {
            mode,
            e: el.e,
            phi0: p.phi0,
            sinx1,
            cosx1,
            akm1,
        })
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
    use super::*;
    use crate::adaptors::transform_xy;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};
    use approx::assert_abs_diff_eq;

    #[test]
    fn proj_stere_stere_ellipsoidal() {
        let p = Proj::from_proj_string("+proj=stere +ellps=GRS80").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222644.85455011716, 110610.88347417387, 0.)),
            ((2., -1., 0.), (222644.85455011716, -110610.88347417528, 0.)),
            ((-2., 1., 0.), (-222644.85455011716, 110610.88347417387, 0.)),
            (
                (-2., -1., 0.),
                (-222644.85455011716, -110610.88347417528, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_stere_stere_spherical() {
        let p = Proj::from_proj_string("+proj=stere +R=6400000").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223407.81025950745, 111737.938996443, 0.)),
            ((2., -1., 0.), (223407.81025950745, -111737.938996443, 0.)),
            ((-2., 1., 0.), (-223407.81025950745, 111737.938996443, 0.)),
            ((-2., -1., 0.), (-223407.81025950745, -111737.938996443, 0.)),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_stere_ups_ellipsoidal() {
        let p = Proj::from_proj_string("+proj=ups +ellps=GRS80").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (2433455.5634384668, -10412543.301512826, 0.)),
            ((2., -1., 0.), (2448749.1185681992, -10850493.419804076, 0.)),
            ((-2., 1., 0.), (1566544.4365615332, -10412543.301512826, 0.)),
            (
                (-2., -1., 0.),
                (1551250.8814318008, -10850493.419804076, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
