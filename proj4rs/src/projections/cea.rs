//!
//! Equal Area Cylindrical Cylindrical (Spherical)
//!
//! Parameters:
//!
//! * lat_ts:
//!
//! See https://proj.org/en/stable/operations/projections/cea.html
//!

use crate::errors::{Error, Result};
use crate::math::consts::{EPS_10, FRAC_PI_2};
use crate::math::{authlat, authset, qsfn};
use crate::parameters::ParamList;
use crate::proj::ProjData;

super::projection! { cea }

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    Sph {
        k0: f64,
    },
    Ell {
        k0: f64,
        e: f64,
        one_es: f64,
        qp: f64,
        apa: (f64, f64, f64),
    },
}

impl Projection {
    pub fn cea(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let (mut k0, t) = match params.try_angular_value("lat_ts")? {
            Some(t) => {
                let k0 = t.cos();
                if k0 < 0. {
                    return Err(Error::InvalidParameterValue(
                        "Invalid value for lat_ts: |lat_ts| should be <= 90\u{00b0}",
                    ));
                }
                (k0, t)
            }
            None => (p.k0, 0.0),
        };

        Ok(if p.ellps.is_ellipsoid() {
            let sint = t.sin();
            k0 /= (1. - p.ellps.es * sint * sint).sqrt();
            Self::Ell {
                k0,
                e: p.ellps.e,
                one_es: p.ellps.one_es,
                qp: qsfn(1., p.ellps.e, p.ellps.one_es),
                apa: authset(p.ellps.es),
            }
        } else {
            Self::Sph { k0 }
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Ell { k0, e, one_es, .. } => {
                Ok((k0 * lam, 0.5 * qsfn(phi.sin(), *e, *one_es) / k0, z))
            }
            Self::Sph { k0 } => Ok((k0 * lam, phi.sin() / k0, z)),
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::Ell { k0, qp, apa, .. } => {
                Ok((x / k0, authlat((2. * y * k0 / qp).asin(), *apa), z))
            }
            Self::Sph { k0 } => {
                let y = y * k0;
                let t = y.abs();
                if t - EPS_10 > 1. {
                    Err(Error::CoordTransOutsideProjectionDomain)
                } else {
                    let phi = if t >= 1. {
                        if y < 0. {
                            -FRAC_PI_2
                        } else {
                            FRAC_PI_2
                        }
                    } else {
                        y.asin()
                    };
                    Ok((x / k0, phi, z))
                }
            }
        }
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

    // lat_ts_0: Lambert cylindrical equal-area

    #[test]
    fn proj_cea_lat_ts_0_e() {
        let p = Proj::from_proj_string("+proj=cea +ellps=GRS80").unwrap();

        // NOTE proj 9 use GRS80 as default ellipsoid

        println!("{:#?}", p.projection());

        let inputs = [(
            (12.09, 47.73, 0.),
            (1345852.643690677360, 4699614.507911851630, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-8);
    }

    #[test]
    fn proj_cea_lat_ts_0_s() {
        let p = Proj::from_proj_string("+proj=cea +R_a +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (12.09, 47.73, 0.),
            (1343596.449131145841, 4711803.232695742510, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-8);
    }

    // lat_ts_30: Berhmann

    #[test]
    fn proj_cea_lat_ts_30_e() {
        let p = Proj::from_proj_string("+proj=cea +lat_ts=30 +ellps=GRS80").unwrap();

        // NOTE proj 9 use GRS80 as default ellipsoid

        println!("{:#?}", p.projection());

        let inputs = [(
            (12.09, 47.73, 0.),
            (1166519.128238123609, 5422104.495923101902, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-8);
    }

    #[test]
    fn proj_cea_lat_ts_30_s() {
        let p = Proj::from_proj_string("+proj=cea +lat_ts=30 +R_a +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (12.09, 47.73, 0.),
            (1163588.657382138772, 5440721.729530871846, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-8);
    }
}
