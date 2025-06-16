//!
//! Miller Cylindrical (Spherical)
//!
//! See https://proj.org/en/stable/operations/projections/mill.html
//!
use crate::errors::Result;
use crate::math::consts::FRAC_PI_4;
use crate::parameters::ParamList;
use crate::proj::ProjData;

super::projection! { mill }

#[derive(Debug, Clone)]
pub(crate) struct Projection {}

impl Projection {
    pub fn mill(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        p.ellps = crate::ellps::Ellipsoid::sphere(p.ellps.a)?;
        Ok(Self {})
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam, (FRAC_PI_4 + phi * 0.4).tan().ln() * 1.25, z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((x, 2.5 * ((0.8 * y).exp().atan() - FRAC_PI_4), z))
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

    #[test]
    fn proj_mill() {
        let p = Proj::from_proj_string("+proj=mill").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (-100., 35.0, 0.),
            (-11131949.079327356070, 4061217.237063715700, 0.),
        )];

        test_proj_forward(&p, &inputs, 1e-8);
        test_proj_inverse(&p, &inputs, 1e-8);
    }
}
