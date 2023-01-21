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
//! lcc, merc, tmerc, utm et aea
//! stere et sterea pour for polar regions.
//!

mod latlong;
mod lcc;

use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::Proj;

pub(crate) type ProjFn = fn(&ProjParams, f64, f64, f64) -> Result<(f64, f64, f64)>;

/// Setup: returned by the init() function
/// Order of members: (params, inverse, forward)
pub(crate) type ProjSetup = (ProjParams, Option<ProjFn>, Option<ProjFn>);

pub(crate) type InitFn = fn(&mut Proj, &ParamList) -> Result<ProjSetup>;

/// Returned by projection looku
pub(crate) struct ProjInit(&'static str, InitFn);

impl ProjInit {
    #[inline(always)]
    pub fn name(&self) -> &'static str {
        self.0
    }

    /// Return a tuple (params, inverse, forward)
    #[inline(always)]
    pub fn init(&self, proj: &mut Proj, params: &ParamList) -> Result<ProjSetup> {
        self.1(proj, params)
    }
}

// Macro for retrieval of parameters from the proj object
// not that makes us writing a match to a unique element.
// XXX Use Into trait instead
macro_rules! downcast {
    ($name:ident, $p:expr) => {
        match $p {
            $crate::projections::ProjParams::$name(data) => data,
            _ => unreachable!(),
        }
    };
}

// Define this in projection definition module
macro_rules! projection {
    ($name:ident) => {
        pub(crate) mod stub {
            use $crate::errors::Result;
            use $crate::parameters::ParamList;
            use $crate::proj::Proj;
            use $crate::projections::{$name, ProjParams, ProjSetup};
            pub(crate) fn init_(p: &mut Proj, params: &ParamList) -> Result<ProjSetup> {
                Ok((
                    ProjParams::$name($name::Projection::init(p, params)?),
                    Some(inverse_),
                    Some(forward_),
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
        }
    };
}

pub(crate) use downcast;
pub(crate) use projection;

macro_rules! declare_projections {
    ($($name:ident),+) => {
        const PROJECTIONS: [ProjInit; 2] = [
        $(
            ProjInit($name::NAME, $name::stub::init_),
        )+
        ];
        #[allow(non_camel_case_types)]
        #[derive(Debug)]
        pub(crate) enum ProjParams {
            NoParams,
            $(
                $name($name::Projection),
            )+
        }
    };
}

// ----------------------------
// Projection list
// ---------------------------

#[rustfmt::skip]
declare_projections! [
    latlong,
    lcc
];

///
/// Return the datum definition
///
pub(crate) fn find_projection(name: &str) -> Option<&ProjInit> {
    PROJECTIONS
        .iter()
        .find(|d| d.name().eq_ignore_ascii_case(name))
}
