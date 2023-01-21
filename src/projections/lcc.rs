//!
//! Lambert Conformal Conic
//!
//! Paramètres:
//!
//! proj: lcc
//!
//! lat_0: the reference latitude
//! lon_0: the reference longitude
//! lat_1: first standard parallel
//! lat_2: second standard parallel
//! x_0: x offset in meters
//! y_0: y offset in meters
//!

use crate::consts::{EPS_10, FRAC_PI_2, FRAC_PI_4};
use crate::errors::{Error, Result};
use crate::math::{msfn, phi2, tsfn};
use crate::parameters::ParamList;
use crate::proj::Proj;

use super::{ProjParams, ProjSetup};

pub(super) const NAME: &str = "lcc";

#[derive(Debug, Default)]
pub(crate) struct Projection {
    phi1: f64,
    phi2: f64,
    n: f64,
    rho0: f64,
    c: f64,
    ellips: bool,
    e: f64,
    k0: f64,
}

impl Projection {
    pub fn init(p: &mut Proj, params: &ParamList) -> Result<ProjSetup> {
        let phi1 = params.try_angular_value("lat_1")?.unwrap_or(0.);
        let phi2 = params.try_angular_value("lat_2")?.unwrap_or_else(|| {
            p.phi0 = p.phi0.or(Some(phi1));
            phi1
        });

        // Standard Parallels cannot be equal and on opposite sides of the equator
        if (phi1 + phi2).abs() < EPS_10 {
            return Err(Error::ProjErrConicLatEqual);
        }

        let phi0 = p.phi0();

        let sinphi = phi1.sin();
        let cosphi = phi1.cos();
        let secant = (phi1 - phi2).abs() >= EPS_10;

        let el = &p.ellps;

        let ellips = el.es != 0.;

        let (mut n, mut c, mut rho0);

        if ellips {
            let m1 = msfn(sinphi, cosphi, el.es);
            let ml1 = tsfn(phi1, sinphi, el.e);
            // secant zone
            n = if secant {
                let sinphi2 = phi2.sin();
                (m1 / msfn(sinphi2, phi2.cos(), el.es)).ln() /
                (ml1 / tsfn(phi2, sinphi2, el.e)).ln()
            } else {
                sinphi
            };
            //rho0 = m1 * ml1.powf(-n);
            //c = rho0 / n;
            c = m1 * ml1.powf(-n) / n;
            rho0 = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
                0.
            } else {
                c * tsfn(phi0, phi0.sin(), el.e).powf(n)
            }
        } else {
            n = if secant {
                (cosphi / phi2.cos()).ln()
                    / ((FRAC_PI_4 + 0.5 * phi2).tan() / (FRAC_PI_4 + 0.5 * phi1).tan()).ln()
            } else {
                sinphi            
            };
            c = cosphi * (FRAC_PI_4 + 0.5 * phi1).tan().powf(n) / n;
            rho0 = if (phi0.abs() - FRAC_PI_2).abs() < EPS_10 {
                0.
            } else {
                c * (FRAC_PI_4 + 0.5 * phi0).tan().powf(-n)
            }
        }

        Ok((
            ProjParams::lcc(Self {
                phi1,
                phi2,
                n,
                rho0,
                c,
                ellips,
                e: el.e,
                k0: p.k0(),
            }),
            Some(Self::inverse),
            Some(Self::forward),
        ))
    }

    pub fn forward(p: &Proj, mut lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let q = super::downcast!(lcc, p);

        let rho = if (phi.abs() - FRAC_PI_2).abs() < EPS_10 {
            if (phi * q.n) <= 0. {
                return Err(Error::ToleranceConditionError);
            } else {
                0.
            }
        } else {
            q.c * if q.ellips {
                tsfn(phi, phi.sin(), q.e).powf(q.n)
            } else {
                (FRAC_PI_4 + 0.5 * phi).tan().powf(-q.n)
            }
        };

        lam *= q.n;

        Ok((
            q.k0 * (rho * lam.sin()),
            q.k0 * (q.rho0 - rho * lam.cos()),
            z,
        ))
    }

    pub fn inverse(p: &Proj, mut x: f64, mut y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let q = super::downcast!(lcc, p);

        x /= q.k0;
        y /= q.k0;

        y = q.rho0 - y;
        // XXX Check this version of hypoth against the
        // one given in proj4
        let mut rho = x.hypot(y);
        let (mut lam, mut phi);
        if rho != 0. {
            if q.n < 0. {
                rho = -rho;
                x = -x;
                y = -y;
            }
            phi = if q.ellips {
                phi2((rho / q.c).powf(1. / q.n), q.e)?
            } else {
                2. * (q.c / rho).powf(1. / q.n).atan() - FRAC_PI_2
            };
            lam = x.atan2(y) / q.n;
        } else {
            lam = 0.;
            phi = if q.n > 0. { FRAC_PI_2 } else { -FRAC_PI_2 };
        }
        Ok((lam, phi, z))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adaptors::transform_xy;
    use crate::proj::Proj;
    use crate::consts::EPS_10;
    use approx::assert_abs_diff_eq;

    fn scale(p: &Proj, xyz: (f64, f64, f64)) -> (f64, f64, f64) {
        (xyz.0 * p.ellps.a, xyz.1 * p.ellps.a, xyz.2)
    }

    fn to_deg(lpz: (f64, f64, f64)) ->(f64, f64, f64) {
        (lpz.0.to_degrees(), lpz.1.to_degrees(), lpz.2)
    }

    fn to_rad(lpz: (f64, f64, f64)) ->(f64, f64, f64) {
        (lpz.0.to_radians(), lpz.1.to_radians(), lpz.2)
    }

    #[test]
    fn proj_lcc_forward() {
        let p = Proj::from_proj_string("+proj=lcc   +ellps=GRS80  +lat_1=0.5 +lat_2=2").unwrap();

        let d = crate::projections::downcast!(lcc, p);

        //println!("{:#?}", d);
        
        let (lam, phi, _) = to_rad((2., 1., 0.));
        let forward = p.forward().unwrap();

        let out = scale(&p, forward(&p, lam, phi, 0.).unwrap());
        assert_eq!(out, (222588.439735968423,  110660.533870799671, 0.));
    }

    #[test]
    fn proj_lcc_inverse() {
        let p = Proj::from_proj_string("+proj=lcc   +ellps=GRS80  +lat_1=0.5 +lat_2=2").unwrap();

        let d = crate::projections::downcast!(lcc, p);

        //println!("{:#?}", d);
    
        let ra = p.ellps.ra;
        // Descale
        let (x, y) = (222588.439735968423 * ra,  110660.533870799671 * ra);
        let inverse = p.inverse().unwrap();

        let out = to_deg(inverse(&p, x, y, 0.).unwrap());
        assert_abs_diff_eq!(out.0, 2., epsilon = EPS_10);
        assert_abs_diff_eq!(out.1, 1., epsilon = EPS_10);
    }

    #[test]
    fn proj_lcc_latlon_to_lcc() {
        let p_from = Proj::from_proj_string("+proj=latlon +ellps=GRS80").unwrap();
        let p_to = Proj::from_proj_string("+proj=lcc +ellps=GRS80 +lat_1=0.5 +lat_2=2").unwrap();

        let (lon_in, lat_in) = (2.0f64.to_radians(), 1.0f64.to_radians());

        let out = transform_xy(&p_from, &p_to, lon_in, lat_in).unwrap();
        assert_eq!(out, (222588.439735968423,  110660.533870799671));
    }

}
