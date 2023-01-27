//!
//! Transverse mercator
//!
//! Provide Evenden/Snyder algorithm
//!
//! Less accurate but slightly faster than the Poder/Engsager algorithm (etmerc)
//!
//! tmerc: "Transverse Mercator" "\n\tCyl, Sph&Ell";
//!
use crate::errors::{Error, Result};
use crate::math::{
    consts::{EPS_10, FRAC_PI_2},
    enfn, inv_mlfn, mlfn, Enfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

#[derive(Debug)]
pub(crate) struct Ell {
    k0: f64,
    es: f64,
    esp: f64,
    ml0: f64,
    en: Enfn,
}

#[derive(Debug)]
pub(crate) struct Sph {
    phi0: f64,
    esp: f64,
    ml0: f64,
}

#[derive(Debug)]
pub(crate) enum Projection {
    Ell(Ell),
    Sph(Sph),
}

use Projection::*;

// Constants for "approximate" transverse mercator
const FC1: f64 = 1.;
const FC2: f64 = 0.5;
const FC3: f64 = 0.16666666666666666666;
const FC4: f64 = 0.08333333333333333333;
const FC5: f64 = 0.05;
const FC6: f64 = 0.03333333333333333333;
const FC7: f64 = 0.02380952380952380952;
const FC8: f64 = 0.01785714285714285714;

impl Projection {
    pub fn estmerc(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        if p.ellps.is_ellipsoid() {
            let es = p.ellps.es;
            let en = enfn(es);
            Ok(Ell(Ell {
                k0: p.k0,
                en,
                es,
                esp: es / (1. - es),
                ml0: mlfn(p.phi0, p.phi0.sin(), p.phi0.cos(), en),
            }))
        } else {
            Ok(Sph(Sph {
                phi0: p.phi0,
                esp: p.k0,
                ml0: 0.5 * p.k0,
            }))
        }
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Ell(e) => e.forward(lam, phi, z),
            Sph(s) => s.forward(lam, phi, z),
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Ell(e) => e.inverse(x, y, z),
            Sph(s) => s.inverse(x, y, z),
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}
// ---------------
// Ellipsoidal
// ---------------
#[rustfmt::skip]
impl Ell {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Fail if our longitude is more than 90 degrees from the
        // central meridian since the results are essentially garbage.
        if lam < -FRAC_PI_2 || lam > FRAC_PI_2 {
            return Err(Error::LatOrLongExceedLimit);
        }

        let (sinphi, cosphi) = phi.sin_cos();
        let mut t = if cosphi.abs() > EPS_10 {
            sinphi / cosphi
        } else {
            0.
        };
        t *= t;
        let mut al = cosphi * lam;
        let als = al * al;
        al /= (1. - self.es * sinphi * sinphi).sqrt();
        let n = self.esp * cosphi * cosphi;
        Ok((
            // x
            self.k0 * al * (FC1 + 
                FC3 * als * (1. - t + n + 
                FC5 * als * (5. + t * (t - 18.) + n * (14. - 58. * t) +
                FC7 * als * (61. + t * (t * (179. - t) - 479.))))),
            // y
            self.k0 * (mlfn(phi, sinphi, cosphi, self.en) - self.ml0 +
                sinphi * al * lam * FC2 * (1. + 
                FC4 * als * (5. - t + n * (9. + 4. * n) + 
                FC6 * als * (61. + t * (t - 58.) + n * (270. - 330. * t) +
                FC8 * als * (1385. + t * (t * (543. - t) - 3111.)))))),
            // z
            z,
        ))
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let phi = inv_mlfn(self.ml0 + y / self.k0, self.es, self.en)?;
        if phi.abs() >= FRAC_PI_2 {
            Ok((0., if y < 0. { -FRAC_PI_2 } else { FRAC_PI_2 }, z))
        } else {
            let (sinphi, cosphi) = phi.sin_cos();
            let mut t = if cosphi.abs() > 1e-10 {
                sinphi / cosphi
            } else {
                0.
            };
            let n = self.esp * cosphi * cosphi;
            let mut con = 1. - self.es * sinphi * sinphi;
            let d = x * con.sqrt() / self.k0;
            con *= t;
            t *= t;
            let ds = d * d;
            Ok((
                // lam
                d * (FC1
                    - ds * FC3 * (1. + 2. * t + n
                    - ds * FC5 * (5. + t * (28. + 24. * t + 8. * n) + 6. * n
                    - ds * FC7 * (61. + t * (662. + t * (1320. + 720. * t))))))
                    / cosphi,
                // phi
                phi - (con * ds / (1. - self.es))
                    * FC2
                    * (1.
                        - ds * FC4 * (5. + t * (3. - 9. * n) + n * (1. - 4. * n)
                        - ds * FC6 * (61. + t * (90. - 252. * n + 45. * t) + 46. * n
                        - ds * FC8 * (1385. + t * (3633. + t * (4095. + 1575. * t)))))),
                // z
                z,
            ))
        }
    }
}

// ---------------
// Spherical
// ---------------
impl Sph {
    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        // Fail if our longitude is more than 90 degrees from the
        // central meridian since the results are essentially garbage.
        if lam < -FRAC_PI_2 || lam > FRAC_PI_2 {
            return Err(Error::LatOrLongExceedLimit);
        }

        let cosphi = phi.cos();
        let mut b = cosphi * lam.sin();
        if (b.abs() - 1.).abs() <= EPS_10 {
            return Err(Error::ToleranceConditionError);
        }

        let x = self.ml0 * ((1. + b) / (1. - b)).ln();
        let mut y = cosphi * lam.cos() / (1. - b * b).sqrt();

        b = y.abs();
        if b >= 1. {
            if (b - 1.) > EPS_10 {
                return Err(Error::ToleranceConditionError);
            } else {
                y = 0.;
            }
        } else {
            y = y.acos();
        }

        if phi < 0. {
            y = -y;
        }
        y = self.esp * (y - self.phi0);

        Ok((x, y, z))
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut h = (x / self.esp).exp();
        let g = 0.5 * (h - 1. / h);

        h = (self.phi0 + y / self.esp).cos();
        let mut phi = ((1. - h * h) / (1. + g * g)).sqrt().asin();

        // Make sure that phi is on the correct hemisphere when false northing is used
        if y < 0. && -phi + self.phi0 < 0.0 {
            phi = -phi
        }

        Ok((if g != 0.0 || h != 0.0 { g.atan2(h) } else { 0. }, phi, z))
    }
}
