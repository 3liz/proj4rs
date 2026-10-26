//! New Zealand Map Grid
//!
//! ref: <https://proj.org/operations/projections/nzmg.html>
//!
//! Fixed to the International 1924 major axis and an origin at 41°S 173°E, so it
//! ignores `lat_0`, `lon_0`, `x_0` and `y_0`.
//!
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { nzmg }

const EPSLN: f64 = 1e-10;
// The literals of proj nzmg.cpp
#[allow(clippy::excessive_precision)]
const SEC5_TO_RAD: f64 = 0.4848136811095359935899141023;
#[allow(clippy::excessive_precision)]
const RAD_TO_SEC5: f64 = 2.062648062470963551564733573;

// Complex coefficients (re, im) of the conformal polynomial
const BF: [(f64, f64); 6] = [
    (0.7557853228, 0.0),
    (0.249204646, 0.003371507),
    (-0.001541739, 0.041058560),
    (-0.10162907, 0.01727609),
    (-0.26623489, -0.36249218),
    (-0.6870983, -1.1651967),
];

const TPSI: [f64; 10] = [
    0.6399175073,
    -0.1358797613,
    0.063294409,
    -0.02526853,
    0.0117879,
    -0.0055161,
    0.0026906,
    -0.001333,
    0.00067,
    -0.00034,
];

const TPHI: [f64; 9] = [
    1.5627014243,
    0.5185406398,
    -0.03333098,
    -0.1052906,
    -0.0368594,
    0.007317,
    0.01220,
    0.00394,
    -0.0013,
];

type Complex = (f64, f64);

#[inline(always)]
fn mul_add(c: Complex, z: Complex, a: Complex) -> Complex {
    (c.0 + z.0 * a.0 - z.1 * a.1, c.1 + z.0 * a.1 + z.1 * a.0)
}

#[inline(always)]
fn mul(z: Complex, a: Complex) -> Complex {
    (z.0 * a.0 - z.1 * a.1, z.0 * a.1 + z.1 * a.0)
}

/// Evaluate the complex polynomial `z * (c[0] + z * (c[1] + ...))`
///
/// See proj `pj_zpoly1`
fn zpoly1(z: Complex, c: &[Complex]) -> Complex {
    let n = c.len() - 1;
    let a = c[..n].iter().rev().fold(c[n], |a, &ci| mul_add(ci, z, a));
    mul(z, a)
}

/// Evaluate the complex polynomial of `zpoly1` and its derivative
///
/// See proj `pj_zpolyd1`
fn zpolyd1(z: Complex, c: &[Complex]) -> (Complex, Complex) {
    let n = c.len() - 1;
    let mut a = c[n];
    let mut b = a;
    for (k, &ci) in c[..n].iter().rev().enumerate() {
        if k > 0 {
            b = mul_add(a, z, b);
        }
        a = mul_add(ci, z, a);
    }
    (mul(z, a), mul_add(a, z, b))
}

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    phi0: f64,
}

impl Projection {
    pub fn nzmg(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        // Force to International major axis
        p.ellps.a = 6378388.0;
        p.ellps.ra = 1. / p.ellps.a;
        p.lam0 = 173f64.to_radians();
        p.phi0 = (-41f64).to_radians();
        p.x0 = 2510000.;
        p.y0 = 6023150.;
        Ok(Self { phi0: p.phi0 })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let phi = (phi - self.phi0) * RAD_TO_SEC5;
        let r = TPSI.iter().rev().fold(0., |r, &c| c + phi * r) * phi;
        let p = zpoly1((r, lam), &BF);
        Ok((p.1, p.0, z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let mut p = (y, x);
        for _ in 0..20 {
            let (f, fp) = zpolyd1(p, &BF);
            let f = (f.0 - y, f.1 - x);
            let den = fp.0 * fp.0 + fp.1 * fp.1;
            let dp = (
                -(f.0 * fp.0 + f.1 * fp.1) / den,
                -(f.1 * fp.0 - f.0 * fp.1) / den,
            );
            p = (p.0 + dp.0, p.1 + dp.1);
            if dp.0.abs() + dp.1.abs() <= EPSLN {
                let phi = TPHI.iter().rev().fold(0., |phi, &c| c + p.0 * phi);
                return Ok((p.1, self.phi0 + p.0 * phi * SEC5_TO_RAD, z));
            }
        }
        Err(Error::ToleranceConditionError)
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

    // Test values from proj builtins.gie
    #[test]
    fn proj_nzmg() {
        let p = Proj::from_proj_string("+proj=nzmg +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (3352675144.747425, -7043205391.100244, 0.)),
            ((2., -1., 0.), (3691989502.779306, -6729069415.332105, 0.)),
            ((-2., 1., 0.), (4099000768.453239, -7863208779.667249, 0.)),
            ((-2., -1., 0.), (4466166927.369976, -7502531736.628605, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-4);

        let inputs = [
            ((175.482086827, -69.422692183, 0.), (200000., 100000., 0.)),
            ((175.756819473, -69.533571088, 0.), (200000., -100000., 0.)),
            ((134.605119233, -61.459995711, 0.), (-200000., 100000., 0.)),
            ((134.333684316, -61.621553676, 0.), (-200000., -100000., 0.)),
        ];
        test_proj_inverse(&p, &inputs, 1e-9);
    }
}
