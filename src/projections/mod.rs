//!
//! Projections
//!
//! We definitely do not want to use trait object, se we use
//! enum for type abstraction.  
//!
//! Instead of writing massive boilerplate for delegation to inner types
//! (there may be dozens of projection !) , we use
//! pointer to associated function. This spare us writing huge  `match`
//! for each fonction call.
//!
//! Most important projections:
//!
//! lcc, merc, tmerc, utm (etmerc) et aea
//! stere et sterea pour for polar regions.
//!

use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::ProjData;

use std::fmt;

pub(crate) type ProjFn = fn(&ProjParams, f64, f64, f64) -> Result<(f64, f64, f64)>;

/// Setup: returned by the init() function
/// Order of members: (params, inverse, forward)
pub(crate) struct ProjDelegate(ProjParams, ProjFn, ProjFn, bool, bool);

impl ProjDelegate {
    #[inline(always)]
    pub fn inverse(&self, u: f64, v: f64, w: f64) -> Result<(f64, f64, f64)> {
        self.1(&self.0, u, v, w)
    }
    #[inline(always)]
    pub fn forward(&self, u: f64, v: f64, w: f64) -> Result<(f64, f64, f64)> {
        self.2(&self.0, u, v, w)
    }

    pub fn has_inverse(&self) -> bool {
        self.3
    }

    pub fn has_forward(&self) -> bool {
        self.4
    }
}

impl fmt::Debug for ProjDelegate {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:#?}", self.0)
    }
}

pub(crate) type InitFn = fn(&mut ProjData, &ParamList) -> Result<ProjDelegate>;

/// Returned by projection looku
pub(crate) struct ProjInit(&'static str, InitFn);

impl ProjInit {
    #[inline(always)]
    pub fn name(&self) -> &'static str {
        self.0
    }

    /// Return a tuple (params, inverse, forward)
    #[inline(always)]
    pub fn init(&self, proj: &mut ProjData, params: &ParamList) -> Result<ProjDelegate> {
        self.1(proj, params)
    }
}

// Macro for retrieval of parameters from the proj object
// not that makes us writing a match to a unique element.
// XXX Use Into trait instead ?
macro_rules! downcast {
    ($name:ident, $p:expr) => {
        match $p {
            $crate::projections::ProjParams::$name(data) => data,
            _ => unreachable!(),
        }
    };
}

//
// Use the following declaration in projection modules
//
// `super::projection!(projection_name);`
//
macro_rules! projection {
    ($name:ident, $id:expr) => {
        pub(crate) mod stub {
            use $crate::errors::Result;
            use $crate::parameters::ParamList;
            use $crate::proj::ProjData;
            use $crate::projections::{$name, ProjDelegate, ProjParams};
            pub(crate) fn init_(p: &mut ProjData, params: &ParamList) -> Result<ProjDelegate> {
                Ok(ProjDelegate(
                    ProjParams::$name($name::Projection::init(p, params)?),
                    inverse_,
                    forward_,
                    $name::Projection::has_inverse(),
                    $name::Projection::has_forward(),
                ))
            }
            pub(crate) fn inverse_(
                p: &ProjParams,
                u: f64,
                v: f64,
                w: f64,
            ) -> Result<(f64, f64, f64)> {
                $crate::projections::downcast!($name, p).inverse(u, v, w)
            }
            pub(crate) fn forward_(
                p: &ProjParams,
                u: f64,
                v: f64,
                w: f64,
            ) -> Result<(f64, f64, f64)> {
                $crate::projections::downcast!($name, p).forward(u, v, w)
            }

            pub const fn name() -> &'static str {
                $id
            }
        }
    };
}

use downcast;
use projection;

macro_rules! declare_projections {
    ($($name:ident),+) => {
        const PROJECTIONS: [ProjInit; 4] = [
        $(
            ProjInit($name::stub::name(), $name::stub::init_),
        )+
        ];
        #[allow(non_camel_case_types)]
        #[derive(Debug)]
        pub(crate) enum ProjParams {
            $(
                $name($name::Projection),
            )+
        }
    };
}

// ----------------------------
// Projection list
// ---------------------------

mod etmerc;
mod latlong;
mod lcc;
mod utm;

#[rustfmt::skip]
declare_projections! [
    latlong,
    lcc,
    etmerc,
    utm
];

///
/// Return the datum definition
///
pub(crate) fn find_projection(name: &str) -> Option<&ProjInit> {
    PROJECTIONS
        .iter()
        .find(|d| d.name().eq_ignore_ascii_case(name))
}
