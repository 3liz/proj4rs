//!
//! Implementation of the aea (Albers Equal Area) projection.
//! and the leac (Lambert Equal Area Conic) projection
//!

// From proj4 PJ_aea.c
//
// Original author:   Gerald Evenden (1995)
//
use crate::errors::{Error, Result};
use crate::math::{
    consts::{EPS_10, EPS_7, FRAC_PI_2},
    msfn, qsfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { aea, leac }

const PHI_NITER: usize = 15;

// determine latitude angle phi1
#[inline]
fn phi1_inv(qs: f64, e: f64, one_es: f64) -> Result<f64> {
    let mut phi = (0.5 * qs).asin();
    if e < EPS_7 {
        Ok(phi)
    } else {
        let mut i = PHI_NITER;
        let (mut sinphi, mut cosphi, mut con, mut com, mut dphi);
        while i > 0 {
            (sinphi, cosphi) = phi.sin_cos();
            con = e * sinphi;
            com = 1. - con * con;
            dphi = 0.5 * com * com / cosphi
                * (qs / one_es - sinphi / com + 0.5 / e * ((1. - con) / (1. + con)).ln());
            phi += dphi;

            if dphi.abs() <= EPS_10 {
                break;
            }

            i -= 1;
        }
        if i == 0 {
            Err(Error::ToleranceConditionError)
        } else {
            Ok(phi)
        }
    }
}

#[derive(Debug, Default, Clone)]
pub(crate) struct Projection {
    e: f64,
    one_es: f64,
    ec: f64,
    n: f64,
    n2: f64,
    c: f64,
    dd: f64,
    rho0: f64,
}

impl Projection {
    pub fn init(p: &ProjData, phi1: f64, phi2: f64) -> Result<Self> {
        if (phi1 + phi2).abs() < EPS_10 {
            return Err(Error::ProjErrConicLatEqual);
        }

        let el = &p.ellps;
        let (sinphi, cosphi) = phi1.sin_cos();
        let mut n = sinphi;
        let secant = (phi1 - phi2).abs() >= EPS_10;

        if el.is_ellipsoid() {
            let m1 = msfn(sinphi, cosphi, el.es);
            let ml1 = qsfn(sinphi, el.e, el.one_es);
            if ml1.is_infinite() {
                return Err(Error::ToleranceConditionError);
            }

            if secant {
                let (sinphi2, cosphi2) = phi2.sin_cos();

                let m2 = msfn(sinphi2, cosphi2, el.es);
                let ml2 = qsfn(sinphi2, el.e, el.one_es);
                if ml2.is_infinite() || ml1 == ml2 {
                    return Err(Error::ToleranceConditionError);
                }
                n = (m1 * m1 - m2 * m2) / (ml2 - ml1);
            }

            let ec = 1. - 0.5 * el.one_es * ((1. - el.e) / (1. + el.e)).ln() / el.e;
            let c = m1 * m1 + n * ml1;
            let dd = 1. / n;
            let n2 = n + n;
            let rho0 = dd * (c - n * qsfn(p.phi0.sin(), el.e, el.one_es)).sqrt();

            Ok(Self {
                e: el.e,
                one_es: el.one_es,
                ec,
                n,
                n2,
                c,
                dd,
                rho0,
            })
        } else {
            if secant {
                n = 0.5 * (n + phi2.sin());
            }
            let dd = 1. / n;
            let n2 = n + n;
            let c = cosphi * cosphi + n2 * sinphi;
            let rho0 = dd * (c - n2 * p.phi0.sin()).sqrt();
            Ok(Self {
                e: el.e,
                one_es: el.one_es,
                ec: 1.,
                n,
                n2,
                c,
                dd,
                rho0,
            })
        }
    }

    #[inline]
    fn is_ellipse(&self) -> bool {
        self.e != 0.
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rho = self.c
            - if self.is_ellipse() {
                self.n * qsfn(phi.sin(), self.e, self.one_es)
            } else {
                self.n2 * phi.sin()
            };

        if rho < 0. {
            Err(Error::ToleranceConditionError)
        } else {
            let rho = self.dd * rho.sqrt();
            let (sin_i, cos_i) = (lam * self.n).sin_cos();
            Ok((rho * sin_i, self.rho0 - rho * cos_i, z))
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (mut xx, mut yy) = (x, self.rho0 - y);
        let mut rho = xx.hypot(yy);
        if rho != 0. {
            if self.n < 0. {
                rho = -rho;
                xx = -xx;
                yy = -yy;
            }
            let mut phi = rho / self.dd;
            if self.is_ellipse() {
                phi = (self.c - phi * phi) / self.n;
                phi = if (self.ec - phi.abs()).abs() > EPS_7 {
                    phi1_inv(phi, self.e, self.one_es)?
                } else if phi < 0. {
                    -FRAC_PI_2
                } else {
                    FRAC_PI_2
                }
            } else {
                phi = (self.c - phi * phi) / self.n2;
                phi = if phi.abs() <= 1. {
                    phi.asin()
                } else if phi < 0. {
                    -FRAC_PI_2
                } else {
                    FRAC_PI_2
                }
            }
            Ok((xx.atan2(yy) / self.n, phi, z))
        } else {
            Ok((0., if self.n > 0. { FRAC_PI_2 } else { -FRAC_PI_2 }, z))
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }

    // ------------
    // aea
    // -----------
    pub fn aea(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::init(
            p,
            params.try_angular_value("lat_1")?.unwrap_or(0.), // phi1
            params.try_angular_value("lat_2")?.unwrap_or(0.), // phi2
        )
    }

    // ----------
    // leac
    // ----------
    pub fn leac(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::init(
            p,
            params.try_angular_value("lat_1")?.unwrap_or(0.), // phi1
            if params.check_option("south")? {
                -FRAC_PI_2
            } else {
                FRAC_PI_2
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_aea_aea_ellipsoidal() {
        let p = Proj::from_proj_string("+proj=aea +ellps=GRS80 +lat_1=0 +lat_2=2").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222571.60875710563, 110653.32674302977, 0.)),
            ((2., -1., 0.), (222706.30650839131, -110484.26714439997, 0.)),
            ((-2., 1., 0.), (-222571.60875710563, 110653.32674302977, 0.)),
            (
                (-2., -1., 0.),
                (-222706.30650839131, -110484.26714439997, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_aea_aea_spherical() {
        let p = Proj::from_proj_string("+proj=aea +R=6400000 +lat_1=0 +lat_2=2").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223334.08517088494, 111780.43188447191, 0.)),
            ((2., -1., 0.), (223470.15499168713, -111610.33943099028, 0.)),
            ((-2., 1., 0.), (-223334.08517088494, 111780.43188447191, 0.)),
            (
                (-2., -1., 0.),
                (-223470.15499168713, -111610.33943099028, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_aea_leac_ellipsoidal() {
        let p = Proj::from_proj_string("+proj=leac +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (220685.14054297868, 112983.50088939646, 0.)),
            ((2., -1., 0.), (224553.31227982609, -108128.63674487274, 0.)),
            ((-2., 1., 0.), (-220685.14054297868, 112983.50088939646, 0.)),
            (
                (-2., -1., 0.),
                (-224553.31227982609, -108128.63674487274, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_aea_leac_spherical() {
        let p = Proj::from_proj_string("+proj=leac +R=6400000").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (221432.86859285168, 114119.45452653214, 0.)),
            ((2., -1., 0.), (225331.72412711097, -109245.82943505641, 0.)),
            ((-2., 1., 0.), (-221432.86859285168, 114119.45452653214, 0.)),
            (
                (-2., -1., 0.),
                (-225331.72412711097, -109245.82943505641, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
