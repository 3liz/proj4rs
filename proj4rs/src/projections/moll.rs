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
#[allow(unused_imports)]
use crate::errors::{Error, Result};
use crate::math::{
    aasin,
    consts::{FRAC_PI_2, PI, TAU},
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { moll, wag4, wag5 }

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    c_x: f64,
    c_y: f64,
    c_p: f64,
    /// p = pi/2: the poles need the special handling in forward() and inverse()
    mollweide: bool,
}

/// s - sin(s), accurate also for small |s|, where the subtraction would lose most of the
/// digits
fn s_minus_sin_s(s: f64) -> f64 {
    if s.abs() < 0.25 {
        let s2 = s * s;
        s * s2
            * (1. / 6.
                - s2 * (1. / 120.
                    - s2 * (1. / 5040.
                        - s2 * (1. / 362880. - s2 * (1. / 39916800. - s2 / 6227020800.)))))
    } else {
        s - s.sin()
    }
}

/// Mollweide close to the poles. With theta = (pi - s) / 2 and chi the colatitude,
/// 2 * theta + sin(2 * theta) = pi * sin(phi) becomes s - sin(s) = 2 * pi * sin(chi / 2)^2,
/// which can be solved accurately where the Newton iteration on theta converges only
/// linearly. Returns s.
fn mollweide_polar_s(chi: f64) -> f64 {
    let sh = (0.5 * chi).sin();
    let w = 2. * PI * sh * sh;
    let mut s = (6. * w).cbrt(); // <= the root, as s - sin(s) <= s^3 / 6
    for _ in 0..20 {
        let sh2 = (0.5 * s).sin();
        let d = 2. * sh2 * sh2; // 1 - cos(s)
        if d == 0. {
            break;
        }
        let ds = (s_minus_sin_s(s) - w) / d;
        s -= ds;
        if ds.abs() <= 1e-12 * s {
            break;
        }
    }
    s
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
            mollweide: false,
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
            mollweide: pp == FRAC_PI_2,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, mut phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        // 30 iterations, as in PROJ (10 were not enough for wag4 near the poles)
        const NITER: isize = 30;
        const TOL: f64 = 1e-7;

        if self.mollweide && phi.abs() > FRAC_PI_2 - 0.5 {
            let s = mollweide_polar_s(FRAC_PI_2 - phi.abs());
            // cos(theta) = sin(s / 2), sin(theta) = cos(s / 2)
            return Ok((
                self.c_x * lam * (0.5 * s).sin(),
                (self.c_y * (0.5 * s).cos()).copysign(phi),
                z,
            ));
        }

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

    /// Close to the poles: cos(s / 2) = |y| / c_y, with s as in mollweide_polar_s(), and the
    /// colatitude chi from 2 * sin(chi / 2)^2 = (s - sin(s)) / pi, instead of asin() of a
    /// value close to 1. Returns None away from the poles.
    #[inline(always)]
    fn polar_inverse(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        if !self.mollweide || y.abs() <= 0.5 * self.c_y || y.abs() > self.c_y {
            return None;
        }
        let half_s = 2. * (0.5 * (self.c_y - y.abs()) / self.c_y).sqrt().asin();
        let den = self.c_x * half_s.sin();
        let lam = if den == 0. { 0. } else { x / den }; // at the pole x is 0
        let chi = 2. * (s_minus_sin_s(2. * half_s) / (2. * PI)).sqrt().asin();
        Some((lam, (FRAC_PI_2 - chi).copysign(y)))
    }

    /// On the map outline the longitude can round to slightly more than pi: accept it
    /// within the given relative tolerance and clamp, so that points on the outline itself
    /// are not rejected
    #[cfg(not(feature = "proj4js-compat"))]
    #[inline(always)]
    fn check_lam(lam: f64, tol: f64) -> Result<f64> {
        if lam.abs() <= PI * (1. + tol) {
            Ok(lam.clamp(-PI, PI))
        } else {
            Err(Error::CoordinateOutOfRange)
        }
    }

    #[cfg(not(feature = "proj4js-compat"))]
    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if let Some((lam, phi)) = self.polar_inverse(x, y) {
            // Close to the poles lam comes from c_y - |y|, whose rounding grows like
            // c_y / (c_y - |y|): the longitude is that ill-conditioned by the input itself
            let tol = 4. * f64::EPSILON * (1. + self.c_y / (self.c_y - y.abs()));
            return Ok((Self::check_lam(lam, tol)?, phi, z));
        }
        let mut phi = aasin(y / self.c_y)?;
        let lam = Self::check_lam(x / (self.c_x * phi.cos()), 4. * f64::EPSILON)?;
        phi += phi;
        phi = aasin((phi + phi.sin()) / self.c_p)?;
        Ok((lam, phi, z))
    }

    /// This is an infaillible version of Mollwey
    /// While proj version is faillible
    #[cfg(feature = "proj4js-compat")]
    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if let Some((lam, phi)) = self.polar_inverse(x, y) {
            return Ok((lam.clamp(-PI, PI), phi, z));
        }
        let mut phi = aasin(y / self.c_y)?;
        let mut lam = x / (self.c_x * phi.cos());
        if lam.abs() > PI {
            lam = PI * lam.signum();
        }
        phi += phi;
        phi = aasin((phi + phi.sin()) / self.c_p)?;
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
