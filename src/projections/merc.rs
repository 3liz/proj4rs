//!
//! Pseudo Mercator
//!
//! merc: "Mercator" "\n\tCyl, Sph&Ell\n\tlat_ts="
//! webmerc: "Web Mercator / Pseudo Mercator" "\n\tCyl, Ell\n\t"
//!

// Projection stub
super::projection! { merc, webmerc }

use crate::errors::{Error, Result};
use crate::math::{
    asinh,
    consts::{EPS_10, FRAC_PI_2},
    msfn, phi2,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

#[derive(Debug)]
pub(crate) struct Projection {
    is_ellps: bool,
    k0: f64,
    e: f64,
}

impl Projection {
    pub fn merc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let phits: Option<f64> = params.try_value("lat_ts")?;
        if let Some(phits) = phits {
            if phits >= FRAC_PI_2 {
                return Err(Error::InvalidParameterValue(
                    "lat_ts larger than 90 degrees",
                ));
            }
        }

        if p.ellps.is_ellipsoid() {
            if let Some(phits) = phits {
                p.k0 = msfn(phits.sin(), phits.cos(), p.ellps.es);
            }
        } else if let Some(phits) = phits {
            p.k0 = phits.cos();
        }

        Ok(Self {
            is_ellps: p.ellps.is_ellipsoid(),
            k0: p.k0,
            e: p.ellps.e,
        })
    }

    pub fn webmerc(p: &mut ProjData, _params: &ParamList) -> Result<Self> {
        p.k0 = 1.0;
        Ok(Self {
            is_ellps: false,
            k0: p.k0,
            e: p.ellps.e,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        if (phi.abs() - FRAC_PI_2).abs() <= EPS_10 {
            return Err(Error::ToleranceConditionError);
        }
        if self.is_ellps {
            let (sphi, cphi) = phi.sin_cos();
            Ok((
                self.k0 * lam,
                self.k0 * (asinh(sphi / cphi) - self.e * (self.e * sphi).atanh()),
                z,
            ))
        } else {
            Ok((self.k0 * lam, self.k0 * asinh(phi.tan()), z))
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.is_ellps {
            Ok((x / self.k0, phi2((-y / self.k0).exp(), self.e)?, z))
        } else {
            Ok((x / self.k0, (y / self.k0).sinh().atan(), z))
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
    use super::*;
    use crate::adaptors::transform_xy;
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};
    use approx::assert_abs_diff_eq;

    #[test]
    fn proj_merc_merc_ellps() {
        let p = Proj::from_proj_string("+proj=merc +ellps=GRS80").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        // XXX: this differs in y from Proj output: 110579.96521824962
        // because of asinh definition: using the same asinh definition
        // leads to same results up to 1e-11m.
        let inputs = [
            ((2., 1., 0.), (222638.98158654713, 110579.96521825077, 0.)),
            ((2., -1., 0.), (222638.98158654713, -110579.96521825077, 0.)),
            ((-2., 1., 0.), (-222638.98158654713, 110579.96521825077, 0.)),
            (
                (-2., -1., 0.),
                (-222638.98158654713, -110579.96521825077, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_merc_merc_sph() {
        let p = Proj::from_proj_string("+proj=merc +R=6400000").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223402.14425527418, 111706.74357494547, 0.)),
            ((2., -1., 0.), (223402.14425527418, -111706.74357494547, 0.)),
            ((-2., 1., 0.), (-223402.14425527418, 111706.74357494547, 0.)),
            (
                (-2., -1., 0.),
                (-223402.14425527418, -111706.74357494547, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_merc_webmerc_ellps() {
        let p = Proj::from_proj_string("+proj=webmerc +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222638.98158654713, 111325.14286638626, 0.)),
            ((2., -1., 0.), (222638.98158654713, -111325.14286638626, 0.)),
            ((-2., 1., 0.), (-222638.98158654713, 111325.14286638626, 0.)),
            (
                (-2., -1., 0.),
                (-222638.98158654713, -111325.14286638626, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_merc_webmerc_sph() {
        let p = Proj::from_proj_string("+proj=webmerc +R=6400000").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223402.14425527418, 111706.74357494547, 0.)),
            ((2., -1., 0.), (223402.14425527418, -111706.74357494547, 0.)),
            ((-2., 1., 0.), (-223402.14425527418, 111706.74357494547, 0.)),
            (
                (-2., -1., 0.),
                (-223402.14425527418, -111706.74357494547, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
