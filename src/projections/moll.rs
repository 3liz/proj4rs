//!
//! Mollweide Pseudocylindrical and derivatives
//!
//! ref: <https://proj.org/operations/projections/moll.html>
//!
//! moll: "Mollweide" "\n\tPCyl., Sph.";
//! wag4: "Wagner IV" "\n\tPCyl., Sph.";
//! wag5: "Wagner V" "\n\tPCyl., Sph.";
//!
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::math::{
    aasin,
    consts::{FRAC_PI_2, PI, TAU},
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { moll, wag4, wag5 }

#[derive(Debug)]
pub(crate) struct Projection {
    c_x: f64,
    c_y: f64,
    c_p: f64,
}

impl Projection {
    pub fn moll(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        Self::new(p, FRAC_PI_2)
    }

    pub fn wag4(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        Self::new(p, PI / 3.)
    }

    pub fn wag5(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        // Map from sphere
        p.ellps = Ellipsoid::sphere(p.ellps.a)?;

        Ok(Self {
            c_x: 0.90977,
            c_y: 1.65014,
            c_p: 3.00896,
        })
    }

    fn new(p: &mut ProjData, pp: f64) -> Result<Self> {
        // Map from sphere
        p.ellps = Ellipsoid::sphere(p.ellps.a)?;

        let p2 = pp + pp;
        let sp = pp.sin();
        let c_p = p2 + p2.sin();
        let r = (TAU * sp / c_p).sqrt();

        Ok(Self {
            c_x: 2. * r / PI,
            c_y: r / sp,
            c_p,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, mut phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        const NITER: isize = 10;
        const TOL: f64 = 1e-7;

        let k = self.c_p * phi.sin();
        let mut i = NITER;
        while i > 0 {
            let v = (phi + phi.sin() - k) / (1. + phi.cos());
            phi -= v;
            if v.abs() < TOL {
                break;
            }
            i -= 1;
        }
        if i == 0 {
            phi = FRAC_PI_2 * phi.signum();
        } else {
            phi *= 0.5;
        }
        Ok((self.c_x * lam * phi.cos(), self.c_y * phi.sin(), z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut phi = aasin(y / self.c_y)?;
        let lam = x / (self.c_x * phi.cos());
        if lam.abs() < PI {
            phi += phi;
            phi = aasin((phi + phi.sin()) / self.c_p)?;
            Ok((lam, phi, z))
        } else {
            Err(Error::CoordinateOutOfRange)
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
    use crate::adaptors::transform_xy;
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_moll() {
        let p = Proj::from_proj_string("+proj=moll").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (200426.67539284358, 123642.46137843542, 0.)),
            ((2., -1., 0.), (200426.67539284358, -123642.46137843542, 0.)),
            ((-2., 1., 0.), (-200426.67539284358, 123642.46137843542, 0.)),
            (
                (-2., -1., 0.),
                (-200426.67539284358, -123642.46137843542, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_wag4() {
        let p = Proj::from_proj_string("+proj=wag4").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (192142.59162431932, 128974.11846682805, 0.)),
            ((2., -1., 0.), (192142.59162431932, -128974.11846682805, 0.)),
            ((-2., 1., 0.), (-192142.59162431932, 128974.11846682805, 0.)),
            (
                (-2., -1., 0.),
                (-192142.59162431932, -128974.11846682805, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_wag5() {
        let p = Proj::from_proj_string("+proj=wag5").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (202532.80926341165, 138177.98447111444, 0.)),
            ((2., -1., 0.), (202532.80926341165, -138177.98447111444, 0.)),
            ((-2., 1., 0.), (-202532.80926341165, 138177.98447111444, 0.)),
            (
                (-2., -1., 0.),
                (-202532.80926341165, -138177.98447111444, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
