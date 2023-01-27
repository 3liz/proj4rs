//!
//! Oblique Stereographic Alternative
//!
//! from PJ_sterea.c (proj 5.2.0)
//!
//! sterea: "Oblique Stereographic Alternative" "\n\tAzimuthal, Sph&Ell"
//!
use crate::errors::{Error, Result};
use crate::math::{gauss, gauss_ini, inv_gauss, Gauss};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { sterea }

#[derive(Debug)]
pub(crate) struct Projection {
    k0: f64,
    phic0: f64,
    cosc0: f64,
    sinc0: f64,
    r2: f64,
    en: Gauss,
}

impl Projection {
    pub fn sterea(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        let (en, phic0, r) = gauss_ini(p.ellps.e, p.phi0)?;
        let (sinc0, cosc0) = phic0.sin_cos();
        Ok(Self {
            k0: p.k0,
            phic0,
            cosc0,
            sinc0,
            r2: 2. * r,
            en,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (lam, phi) = gauss(lam, phi, &self.en);
        let (sinc, cosc) = phi.sin_cos();
        let cosl = lam.cos();
        let k = self.k0 * self.r2 / (1. + self.sinc0 * sinc + self.cosc0 * cosc * cosl);
        Ok((
            k * cosc * lam.sin(),
            k * (self.cosc0 * sinc - self.sinc0 * cosc * cosl),
            z,
        ))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let x = x / self.k0;
        let y = y / self.k0;
        let rho = x.hypot(y);
        let (lam, phi) = if rho != 0.0 {
            let c = 2. * rho.atan2(self.r2);
            let (sinc, cosc) = c.sin_cos();
            inv_gauss(
                (x * sinc).atan2(rho * self.cosc0 * cosc - y * self.sinc0 * sinc),
                (cosc * self.sinc0 + y * sinc * self.cosc0 / rho).asin(),
            &self.en)
        } else {
            inv_gauss(0., self.phic0, &self.en)
        }?;
        Ok((lam, phi, z))
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
    fn proj_sterea() {
        let p = Proj::from_proj_string("+proj=sterea +ellps=GRS80").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222644.89410919772, 110611.09187173686, 0.)),
            ((2., -1., 0.), (222644.89410919772, -110611.09187173827, 0.)),
            ((-2., 1., 0.), (-222644.89410919772, 110611.09187173686, 0.)),
            (
                (-2., -1., 0.),
                (-222644.89410919772, -110611.09187173827, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
