//!
//! Derived values for ellipsoids
//!
//! A complete ellipsoid definition comprises a size (primary) and a shape (secondary)
//! parameter.
//!
//! Size parameters supported are:
//!     R, defining the radius of a spherical planet
//!     a, defining the semimajor axis of an ellipsoidal planet
//!
//! Shape parameters supported are:
//!     rf, the reverse flattening of the ellipsoid
//!     f,  the flattening of the ellipsoid
//!     es, the eccentricity squared
//!     e,  the eccentricity
//!     b,  the semiminor axis
//!
//! The ellps=xxx parameter provides both size and shape for a number of built in
//! ellipsoid definitions.
//!
//! The ellipsoid definition may be augmented with a spherification flag, turning
//! the ellipsoid into a sphere with features defined by the ellipsoid.
//!
//! Spherification parameters supported are:
//!     R_A, which gives a sphere with the same surface area as the ellipsoid
//!     R_V, which gives a sphere with the same volume as the ellipsoid
//!     R_a, which gives a sphere with R = (a + b)/2   (arithmetic mean)
//!     R_g, which gives a sphere with R = sqrt(a*b)   (geometric mean)
//!     R_h, which gives a sphere with R = 2*a*b/(a+b) (harmonic mean)
//!     R_lat_a=phi, which gives a sphere with R being the arithmetic mean of
//!         of the corresponding ellipsoid at latitude phi.
//!     R_lat_g=phi, which gives a sphere with R being the geometric mean of
//!         of the corresponding ellipsoid at latitude phi.
//!
#![allow(non_upper_case_globals)]

use crate::ellipsoids::{EllipsoidDefn, FlatteningParam};
use crate::errors::{Error, Result};
use crate::math::consts::EPS_10;
use crate::parameters::ParamList;

use std::ops::ControlFlow;

// series coefficients for calculating ellipsoid-equivalent spheres
const SIXTH: f64 = 1. / 6.;
const RA4: f64 = 17. / 360.;
const RA6: f64 = 67. / 3024.;
const RV4: f64 = 5. / 72.;
const RV6: f64 = 55. / 1296.;

const TOK_rf: &str = "rf";
const TOK_f: &str = "f";
const TOK_es: &str = "es";
const TOK_e: &str = "e";
const TOK_b: &str = "b";

const TOK_R_A: &str = "R_A";
const TOK_R_V: &str = "R_V";
const TOK_R_a: &str = "R_a";
const TOK_R_g: &str = "R_g";
const TOK_R_h: &str = "R_h";

/// A shape parameter
#[allow(non_camel_case_types)]
enum Shape {
    SP_rf(f64),
    SP_f(f64),
    SP_es(f64),
    SP_e(f64),
    SP_b(f64),
}

#[derive(Clone, Debug)]
pub struct Ellipsoid {
    // The linear parameters
    pub a: f64,  // semimajor axis (radius if eccentricity==0)
    pub b: f64,  // semiminor axis
    pub ra: f64, // 1/a
    pub rb: f64, // 1/b

    // The eccentricities
    //pub alpha: f64,   // angular eccentricity
    pub e: f64,  // first  eccentricity
    pub es: f64, // first  eccentricity squared
    //pub e2: f64,    // second eccentricity
    //pub e2s: f64,   // second eccentricity squared
    //pub e3: f64,    // third  eccentricity
    //pub e3s: f64,   // third  eccentricity squared
    pub one_es: f64,  // 1 - e^2
    pub rone_es: f64, // 1/one_es

    // The flattenings
    pub f: f64, // first  flattening
    //pub f2: f64,  // second flattening
    //pub n: f64,   // third  flattening
    pub rf: f64, // 1/f

                 /*
                     pub rf2: f64, // 1/f2
                     pub rn: f64,  // 1/n

                     // This one's for GRS80
                     pub jform: f64, // Dynamic form factor

                     pub es_orig: f64, // es and a before any +proj related adjustment
                     pub a_orig: f64,
                 */
}

use Shape::*;

impl Ellipsoid {
    #[inline]
    pub fn is_sphere(&self) -> bool {
        self.es == 0.
    }

    #[inline]
    pub fn is_ellipsoid(&self) -> bool {
        self.es != 0.
    }

    /// Create sphere
    pub fn sphere(radius: f64) -> Result<Self> {
        if !(radius.is_normal() && radius > 0.) {
            return Err(Error::InvalidParameterValue("Invalid radius"));
        }
        Ok(Self {
            a: radius,
            b: radius,
            ra: 1. / radius,
            rb: 1. / radius,
            e: 0.,
            es: 0.,
            f: 0.,
            rf: f64::INFINITY,
            one_es: 1.,
            rone_es: 1.,
        })
    }

    #[cfg(test)]
    pub fn try_from_ellipsoid(defn: &EllipsoidDefn) -> Result<Self> {
        Self::calc_ellipsoid_params(
            defn.a,
            match defn.rf_or_b {
                FlatteningParam::MinorAxis(b) => SP_b(b),
                FlatteningParam::InvFlat(rf) => SP_rf(rf),
            },
        )
    }

    /// Create ellipsoid from definition and parameters
    pub fn try_from_ellipsoid_with_params(
        defn: &EllipsoidDefn,
        params: &ParamList,
    ) -> Result<Self> {
        // Override "a" ?
        let a = if let Some(p) = params.get("a") {
            p.try_into()?
        } else {
            defn.a
        };
        // Get the shape parameter
        let sp = Self::find_shape_parameter(params).unwrap_or(Ok(match defn.rf_or_b {
            FlatteningParam::MinorAxis(b) => SP_b(b),
            FlatteningParam::InvFlat(rf) => SP_rf(rf),
        }))?;
        Self::calc_ellipsoid_params(a, sp).and_then(|ellps| ellps.spherification(params))
    }

    fn find_shape_parameter(params: &ParamList) -> Option<Result<Shape>> {
        // Shape parameters tokens in order of precedence
        const SHAPE_TOKENS: &[&str] = &[TOK_rf, TOK_f, TOK_es, TOK_e, TOK_b];
        SHAPE_TOKENS.iter().find_map(|tok| {
            params.get(tok).map(|p| {
                p.try_into().map(|v| match *tok {
                    TOK_rf => SP_rf(v),
                    TOK_f => SP_f(v),
                    TOK_es => SP_es(v),
                    TOK_e => SP_e(v),
                    TOK_b => SP_b(v),
                    _ => unreachable!(),
                })
            })
        })
    }

    /// Calculate parameters and return a new ellipsoid
    /// This is the true constructor
    fn calc_ellipsoid_params(a: f64, sp: Shape) -> Result<Self> {
        if a <= 0. {
            return Err(Error::InvalidParameterValue("Invalid major axis"));
        }

        let (mut f, mut rf, mut es, mut e, mut b);
        // We could have return directly a tuple from the match expression
        // but that makes the code less readable and the compiler will check
        // uninitialized variables anyway.
        match sp {
            SP_rf(p_rf) => {
                if p_rf <= 1. {
                    return Err(Error::InvalidParameterValue(
                        "Inverse flattening lower than 1.",
                    ));
                }
                rf = p_rf;
                f = 1. / rf;
                es = 2. * f - f * f;
                e = es.sqrt();
                b = (1.0 - f) * a;
            }
            SP_f(p_f) => {
                if !(0. ..1.).contains(&p_f) {
                    return Err(Error::InvalidParameterValue("Flattening not in [0..1["));
                }
                f = p_f;
                es = 2. * f - f * f;
                e = es.sqrt();
                b = (1.0 - f) * a;
                rf = if f > 0. { 1. / f } else { f64::INFINITY }
            }
            SP_es(p_es) => {
                if !(0. ..1.).contains(&p_es) {
                    return Err(Error::InvalidParameterValue(
                        "Square eccentricity not in [0..1[",
                    ));
                }
                es = p_es;
                e = es.sqrt();
                f = 1. - e.asin().cos();
                b = (1.0 - f) * a;
                rf = if f > 0. { 1. / f } else { f64::INFINITY }
            }
            SP_e(p_e) => {
                if !(0. ..1.).contains(&p_e) {
                    return Err(Error::InvalidParameterValue("Eccentricity not in [0..1["));
                }
                e = p_e;
                es = e * e;
                f = 1. - e.asin().cos();
                b = (1.0 - f) * a;
                rf = if f > 0. { 1. / f } else { f64::INFINITY }
            }
            SP_b(p_b) => {
                if !(p_b > 0. && p_b <= a) {
                    return Err(Error::InvalidParameterValue("Invalid minor axis"));
                }
                b = p_b;
                let a2 = a * a;
                let b2 = b * b;
                es = (a2 - b2) / a2;
                e = es.sqrt();
                f = (a - b) / b;
                rf = if f > 0. { 1. / f } else { f64::INFINITY }
            }
        }

        if (a - b).abs() < EPS_10 {
            b = a;
            es = 0.;
            e = 0.;
            f = 0.;
            rf = f64::INFINITY;
        }

        let one_es = 1. - es;

        Ok(Self {
            a,
            b,
            ra: 1. / a,
            rb: 1. / b,
            e,
            es,
            f,
            rf,
            one_es,
            rone_es: 1. / one_es,
        })
    }

    fn spherification(self, params: &ParamList) -> Result<Self> {
        // Spherification parameter
        const SPHERE_TOKENS: &[&str] = &[TOK_R_A, TOK_R_V, TOK_R_a, TOK_R_g, TOK_R_h];
        match SPHERE_TOKENS.iter().try_for_each(|tok| {
            if params.get(tok).is_some() {
                let es = self.es;
                let a = match *tok {
                    // a sphere with same area as ellipsoid
                    TOK_R_A => 1. - es * (SIXTH + es * (RA4 + es * RA6)),
                    // a sphere with same volume as ellipsoid
                    TOK_R_V => 1. - es * (SIXTH + es * (RV4 + es * RV6)),
                    // a sphere with R = the arithmetic mean of the ellipsoid
                    TOK_R_a => (self.a + self.b) / 2.,
                    // a sphere with R = the geometric mean of the ellipsoid
                    TOK_R_g => (self.a + self.b).sqrt(),
                    // a sphere with R = the harmonic mean of the ellipsoid
                    TOK_R_h => (2. * self.a * self.b) / (self.a + self.b),
                    _ => unreachable!(),
                };
                // Update ellipsoid parameters
                ControlFlow::Break(Self::calc_ellipsoid_params(a, SP_es(0.)))
            } else {
                ControlFlow::Continue(())
            }
        }) {
            ControlFlow::Break(rv) => rv,
            _ => Ok(self),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ellipsoids::constants::*;
    use crate::projstring;

    #[test]
    fn ellps_from_defn() {
        let ellps = Ellipsoid::try_from_ellipsoid(&WGS84).unwrap();

        assert_eq!(ellps.a, 6_378_137.);
        assert_eq!(ellps.rf, 298.257_223_563);
    }

    #[test]
    fn ellps_from_defn_and_params() {
        let ellps = Ellipsoid::try_from_ellipsoid_with_params(
            &WGS84,
            &projstring::parse("+a=6370997.").unwrap(),
        )
        .unwrap();

        assert_eq!(ellps.a, 6_370_997.);
        assert_eq!(ellps.rf, 298.257_223_563);
    }

    fn assert_sphere(ellps: Ellipsoid) {
        assert_eq!(ellps.a, ellps.b);
        assert_eq!(ellps.f, 0.);
        assert_eq!(ellps.rf, f64::INFINITY);
        assert_eq!(ellps.e, 0.);
        assert_eq!(ellps.es, 0.);
    }

    #[test]
    fn ellps_sphere() {
        let ellps = Ellipsoid::sphere(6_378_388.).unwrap();

        assert_eq!(ellps.a, 6_378_388.);
        assert_sphere(ellps);
    }

    #[test]
    fn ellps_from_defn_and_es_zero() {
        let ellps = Ellipsoid::try_from_ellipsoid_with_params(
            &WGS84,
            &projstring::parse("+es=0.").unwrap(),
        )
        .unwrap();

        assert_sphere(ellps);
    }

    #[test]
    fn ellps_spherification() {
        let ellps =
            Ellipsoid::try_from_ellipsoid_with_params(&WGS84, &projstring::parse("+R_A").unwrap())
                .unwrap();

        assert_sphere(ellps);
    }

    #[test]
    fn ellps_invalid_params() {
        fn from_projstring(s: &str) -> Result<Ellipsoid> {
            Ellipsoid::try_from_ellipsoid_with_params(&WGS84, &projstring::parse(s).unwrap())
        }

        assert!(from_projstring("+a=-0.").is_err());
        assert!(from_projstring("+a=-2.").is_err());
        assert!(from_projstring("+es=-1.").is_err());
        assert!(from_projstring("+f=20.").is_err());
    }
}
