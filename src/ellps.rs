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
//!
//!     R_a, which gives a sphere with R = (a + b)/2   (arithmetic mean)
//!     R_g, which gives a sphere with R = sqrt(a*b)   (geometric mean)
//!     R_h, which gives a sphere with R = 2*a*b/(a+b) (harmonic mean)
//!
//!     R_lat_a=phi, which gives a sphere with R being the arithmetic mean of
//!         of the corresponding ellipsoid at latitude phi.
//!     R_lat_g=phi, which gives a sphere with R being the geometric mean of
//!         of the corresponding ellipsoid at latitude phi.
//!

use crate::constants::EPSLN;
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::{datums, ellipsoids, projstring};

// series coefficients for calculating ellipsoid-equivalent spheres
const SIXTH: f64 = 1. / 6.;
const RA4: f64 = 17. / 360.;
const RA6: f64 = 67. / 3024.;
const RV4: f64 = 5. / 72.;
const RV6: f64 = 55. / 1296.;

#[derive(Default, Clone)]
pub(crate) struct PJConsts {
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
                 //pub rf2: f64, // 1/f2
                 //pub rn: f64,  // 1/n

                 // This one's for GRS80
                 // pub jform: f64, // Dynamic form factor

                 //pub es_orig: f64, // es and a before any +proj related adjustment
                 //pub a_orig: f64,
}

/// A shape parameter
/// by order of precedence
#[allow(non_camel_case_types)]
enum ShapeParameter {
    SP_rf(f64),
    SP_f(f64),
    SP_es(f64),
    SP_e(f64),
    SP_b(f64),
}

use ShapeParameter::*;

impl PJConsts {
    fn _sphere(mut self, radius: f64) -> Self {
        assert!(radius.is_normal());
        assert!(radius > 0.);
        self.a = radius;
        self.b = self.a;
        self.rf = f64::INFINITY;
        self
    }

    /// Crate sphere parameters
    pub fn sphere(radius: f64) -> Self {
        Self::default()._sphere(radius)
    }

    /// Calculate parameters given a and es
    ///
    /// Precedence of shape parameters are
    /// "rf", "f", "es", "e", "b"
    fn calc_ellipsoid_params(&mut self, sp: ShapeParameter) -> Result<()> {
        if self.a <= 0. {
            return Err(Error::InvalidParameterValue("Invvalid major axis"));
        }

        let a = self.a;

        match sp {
            SP_rf(rf) => {
                if !(rf >= 0. && rf > 1.) {
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
                if !(f >= 0. && f < 1.) {
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
                if !(es >= 0. && es < 1.) {
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
                if !(e >= 0. && e < 1.) {
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

        if (a - self.b).abs() < EPSLN {
            self.b = a;
            self.es = 0.;
            self.e = 0.;
            self.f = 0.;
            self.rf = f64::INFINITY;
        }

        self.ra = 1. / self.a;
        self.rb = 1. / self.b;

        Ok(())
    }

    /*
            if self.r_a {
                self.a *= 1 - es * (SIXTH + es * (RA4 + es * RA6));
                a2 = self.a * self.a;
                es = 0;
            }
    */
}
