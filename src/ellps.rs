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

use crate::constants::EPSLN;
use crate::ellipsoids::{EllipsoidDefn, FlatteningParam};
use crate::errors::{Error, Result};
use crate::parameters::ParamList;

use std::ops::ControlFlow;

// series coefficients for calculating ellipsoid-equivalent spheres
const SIXTH: f64 = 1. / 6.;
const RA4: f64 = 17. / 360.;
const RA6: f64 = 67. / 3024.;
const RV4: f64 = 5. / 72.;
const RV6: f64 = 55. / 1296.;

#[derive(Default, Clone, Debug)]
pub(crate) struct Ellipsoid {
    // The linear parameters
    pub a: f64,  // semimajor axis (radius if eccentricity==0)
    pub b: f64,  // semiminor axis
    pub ra: f64, // 1/a
    pub rb: f64, // 1/b

    // The eccentricities
    //pub alpha: f64,   // angular eccentricity
    pub e: f64,  // first  eccentricity
    pub es: f64, // first  eccentricity squared
    //pub e2: f64,      // second eccentricity
    //pub e2s: f64,     // second eccentricity squared
    //pub e3: f64,      // third  eccentricity
    //pub e3s: f64,     // third  eccentricity squared
    //pub one_es: f64,  // 1 - e^2
    //pub rone_es: f64, // 1/one_es

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
    pub ep2: f64,
}

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

use Shape::*;

impl Ellipsoid {
    /// Create sphere parameters
    #[inline]
    pub fn sphere(radius: f64) -> Result<Self> {
        Self::default()._sphere(radius)
    }

    fn _sphere(mut self, radius: f64) -> Result<Self> {
        if !(radius.is_normal() && radius > 0.) {
            return Err(Error::InvalidParameterValue("Invalid radius"));
        }
        self.a = radius;
        self.b = self.a;
        self.rf = f64::INFINITY;
        Ok(self)
    }

    /// Create ellipsoid from definition
    #[inline]
    pub fn from_ellipsoid(defn: &EllipsoidDefn) -> Result<Self> {
        Self::default()._from_ellipsoid(defn)
    }

    fn _from_ellipsoid(mut self, defn: &EllipsoidDefn) -> Result<Self> {
        self.a = defn.a;
        self.calc_ellipsoid_params(match defn.rf_or_b {
            FlatteningParam::MinorAxis(b) => SP_b(b),
            FlatteningParam::InvFlat(rf) => SP_rf(rf),
        })?;
        Ok(self)
    }

    /// Create ellipsoid from definition and parameters
    #[inline]
    pub fn from_ellipsoid_with_params(defn: &EllipsoidDefn, params: &ParamList) -> Result<Self> {
        Self::default()._ellipsoid_with_params(defn, params)
    }

    fn _ellipsoid_with_params(mut self, defn: &EllipsoidDefn, params: &ParamList) -> Result<Self> {
        self.a = defn.a;
        // Override "a"
        if let Some(p) = params.get("a") {
            self.a = p.try_into()?;
        }
        // Get the shape parameter
        let sp = Self::find_shape_parameter(params).unwrap_or(Ok(match defn.rf_or_b {
            FlatteningParam::MinorAxis(b) => SP_b(b),
            FlatteningParam::InvFlat(rf) => SP_rf(rf),
        }))?;
        self.calc_ellipsoid_params(sp)?;
        self._spherification(params)?;
        Ok(self)
    }

    fn find_shape_parameter(params: &ParamList) -> Option<Result<Shape>> {
        // Shape parameters tokens in order of precedence
        const SHAPE_TOKENS: &[&str] = &[TOK_rf, TOK_f, TOK_es, TOK_e, TOK_b];
        SHAPE_TOKENS.iter().find_map(|tok| {
            if let Some(p) = params.get(tok) {
                Some(p.try_convert::<f64>().map(|v| match *tok {
                    TOK_rf => SP_rf(v),
                    TOK_f => SP_f(v),
                    TOK_es => SP_es(v),
                    TOK_e => SP_e(v),
                    TOK_b => SP_b(v),
                    _ => unreachable!(),
                }))
            } else {
                None
            }
        })
    }

    /// Calculate parameters
    fn calc_ellipsoid_params(&mut self, sp: Shape) -> Result<()> {
        if self.a <= 0. {
            return Err(Error::InvalidParameterValue("Invalid major axis"));
        }

        let a = self.a;

        match sp {
            SP_rf(rf) => {
                if rf <= 1. {
                    return Err(Error::InvalidParameterValue("Invalid inverse flattening"));
                }
                let f = 1. / rf;
                self.f = f;
                self.rf = rf;
                self.es = 2. * f - f * f;
                self.e = self.es.sqrt();
                self.b = (1.0 - f) * a;
            }
            SP_f(f) => {
                if !(0. ..1.).contains(&f) {
                    return Err(Error::InvalidParameterValue("Invalid flattening"));
                }
                self.f = f;
                self.es = 2. * f - f * f;
                self.e = self.es.sqrt();
                self.b = (1.0 - f) * a;
                if f > 0. {
                    self.rf = 1. / f;
                }
            }
            SP_es(es) => {
                if !(0. ..1.).contains(&es) {
                    return Err(Error::InvalidParameterValue("Invalid eccentricity"));
                }
                self.es = es;
                self.e = es.sqrt();
                self.f = 1. - self.e.asin().cos();
                self.b = (1.0 - self.f) * a;
                if self.f > 0. {
                    self.rf = 1. / self.f;
                }
            }
            SP_e(e) => {
                if !(0. ..1.).contains(&e) {
                    return Err(Error::InvalidParameterValue("Invalid eccentricity"));
                }
                self.es = e * e;
                self.e = e;
                self.f = 1. - self.e.asin().cos();
                self.b = (1.0 - self.f) * a;
                if self.f > 0. {
                    self.rf = 1. / self.f;
                }
            }
            SP_b(b) => {
                if !(b >= 0. && b < a) {
                    return Err(Error::InvalidParameterValue("Invalid minor axis"));
                }
                let a2 = a * a;
                let b2 = b * b;
                self.b = b;
                self.es = (a2 - b2) / a2;
                self.e = self.es.sqrt();
                self.f = (a - b) / b;
                if self.f > 0. {
                    self.rf = 1. / self.f;
                }
            }
        }

        let b = self.b;
        if (a - b).abs() < EPSLN {
            self.b = a;
            self.es = 0.;
            self.e = 0.;
            self.f = 0.;
            self.rf = f64::INFINITY;
            self.ep2 = 0.;
        } else {
            let b2 = b * b;
            self.ep2 = ((a * a) - b2) / b2;
        }

        self.ra = 1. / a;
        self.rb = 1. / b;

        Ok(())
    }

    fn _spherification(&mut self, params: &ParamList) -> Result<()> {
        // Spherification parameter
        const SPHERE_TOKENS: &[&str] = &[TOK_R_A, TOK_R_V, TOK_R_a, TOK_R_g, TOK_R_h];
        match SPHERE_TOKENS.iter().try_for_each(|tok| {
            if params.get(tok).is_some() {
                let es = self.es;
                self.a = match *tok {
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
                ControlFlow::Break(self.calc_ellipsoid_params(SP_es(0.)))
            } else {
                ControlFlow::Continue(())
            }
        }) {
            ControlFlow::Break(rv) => rv,
            _ => Ok(()),
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
        let ellps = Ellipsoid::from_ellipsoid(&WGS84).unwrap();

        assert_eq!(ellps.a, 6_378_137.);
        assert_eq!(ellps.rf, 298.257_223_563);
    }

    #[test]
    fn ellps_from_defn_and_params() {
        let ellps = Ellipsoid::from_ellipsoid_with_params(
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
        let ellps =
            Ellipsoid::from_ellipsoid_with_params(&WGS84, &projstring::parse("+es=0.").unwrap())
                .unwrap();

        assert_sphere(ellps);
    }

    #[test]
    fn ellps_spherification() {
        let ellps =
            Ellipsoid::from_ellipsoid_with_params(&WGS84, &projstring::parse("+R_A").unwrap())
                .unwrap();

        assert_sphere(ellps);
    }

    #[test]
    fn ellps_invalid_params() {
        fn from_projstring(s: &str) -> Result<Ellipsoid> {
            Ellipsoid::from_ellipsoid_with_params(&WGS84, &projstring::parse(s).unwrap())
        }

        assert!(from_projstring("+a=-0.").is_err());
        assert!(from_projstring("+a=-2.").is_err());
        assert!(from_projstring("+es=-1.").is_err());
        assert!(from_projstring("+f=20.").is_err());
    }
}
