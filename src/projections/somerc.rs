//!
//! Swiss Oblique Mercator
//!
//! ref: <https://proj.org/operations/projections/somerc.html>
//!
//! somerc: "Swiss. Obl. Mercator" "\n\tCyl, Ell\n\tFor CH1903";
//!
use crate::errors::{Error, Result};
use crate::math::{
    aasin,
    consts::{EPS_10, FRAC_PI_2, FRAC_PI_4},
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { somerc }

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    e: f64,
    rone_es: f64,
    k: f64,
    c: f64,
    hlf_e: f64,
    k_r: f64,
    cosp0: f64,
    sinp0: f64,
}

#[allow(non_snake_case)]
impl Projection {
    pub fn somerc(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        let el = &p.ellps;
        let hlf_e = 0.5 * el.e;

        let (sinphi, cosphi) = p.phi0.sin_cos();

        let cp = cosphi * cosphi;
        let c = (1. + el.es * cp * cp * el.rone_es).sqrt();
        let sinp0 = sinphi / c;
        let phip0 = aasin(sinp0)?;
        let cosp0 = phip0.cos();
        let sp = sinphi * el.e;
        let k = (FRAC_PI_4 + 0.5 * phip0).tan().ln()
            - c * ((FRAC_PI_4 + 0.5 * p.phi0).tan().ln() - hlf_e * ((1. + sp) / (1. - sp)).ln());
        let k_r = p.k0 * el.one_es.sqrt() / (1. - sp * sp);
        Ok(Self {
            e: el.e,
            rone_es: el.rone_es,
            k,
            c,
            hlf_e,
            k_r,
            cosp0,
            sinp0,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let sp = self.e * phi.sin();
        let phip = 2.
            * ((self.c
                * ((FRAC_PI_4 + 0.5 * phi).tan().ln()
                    - self.hlf_e * ((1. + sp) / (1. - sp)).ln())
                + self.k)
                .exp())
            .atan()
            - FRAC_PI_2;

        let lamp = self.c * lam;
        let cp = phip.cos();
        let phipp = aasin(self.cosp0 * phip.sin() - self.sinp0 * cp * lamp.cos())?;
        let lampp = aasin(cp * lamp.sin() / phipp.cos())?;

        Ok((
            self.k_r * lampp,
            self.k_r * (FRAC_PI_4 + 0.5 * phipp).tan().ln(),
            z,
        ))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        const NITER: isize = 6;

        let phipp = 2. * (((y / self.k_r).exp()).atan() - FRAC_PI_4);
        let lampp = x / self.k_r;
        let cp = phipp.cos();
        let mut phip = aasin(self.cosp0 * phipp.sin() + self.sinp0 * cp * lampp.cos())?;
        let lamp = aasin(cp * lampp.sin() / phip.cos())?;
        let con = (self.k - (FRAC_PI_4 + 0.5 * phip).tan().ln()) / self.c;

        let mut i = NITER;
        while i > 0 {
            let esp = self.e * phip.sin();
            let delp = (con + (FRAC_PI_4 + 0.5 * phip).tan().ln()
                - self.hlf_e * ((1. + esp) / (1. - esp)).ln())
                * (1. - esp * esp)
                * phip.cos()
                * self.rone_es;
            phip -= delp;
            if delp.abs() < EPS_10 {
                break;
            }
            i -= 1;
        }
        if i <= 0 {
            Err(Error::ToleranceConditionError)
        } else {
            Ok((lamp / self.c, phip, z))
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
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_somerc_el() {
        let p = Proj::from_proj_string("+proj=somerc +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222638.98158654713, 110579.96521824898, 0.)),
            ((2., -1., 0.), (222638.98158654713, -110579.96521825089, 0.)),
            ((-2., 1., 0.), (-222638.98158654713, 110579.96521824898, 0.)),
            (
                (-2., -1., 0.),
                (-222638.98158654713, -110579.96521825089, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_somerc_sp() {
        let p = Proj::from_proj_string("+proj=somerc +a=6400000").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223402.14425527418, 111706.74357494408, 0.)),
            ((2., -1., 0.), (223402.14425527418, -111706.74357494518, 0.)),
            ((-2., 1., 0.), (-223402.14425527418, 111706.74357494408, 0.)),
            (
                (-2., -1., 0.),
                (-223402.14425527418, -111706.74357494518, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
