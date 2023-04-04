//!
//! Implemented projections
//!

// We definitely do not want to use trait object, se we use
// enum for type abstraction.
//
// Instead of writing massive boilerplate for delegation to inner types
// (there may be dozens of projection !) , we use
// pointer to associated function. This spare us writing huge  `match`
// for each fonction call.
//
// Most important projections:
//
// lcc, merc, tmerc, utm (etmerc) et aea
// stere et sterea pour for polar regions.
//

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

/// Returned by projection lookup
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
macro_rules! projection_delegate {
    ($name:ident, $($init:ident),+ $(,)?) => {
        pub(crate) mod stub {
            use $crate::errors::Result;
            use $crate::parameters::ParamList;
            use $crate::proj::ProjData;
            use $crate::projections::{$name, ProjDelegate, ProjParams};
            $(pub(crate) fn $init(p: &mut ProjData, params: &ParamList) -> Result<ProjDelegate> {
                Ok(ProjDelegate(
                    ProjParams::$name($name::Projection::$init(p, params)?),
                    inverse_,
                    forward_,
                    $name::Projection::has_inverse(),
                    $name::Projection::has_forward(),
                ))
            })+
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

macro_rules! projection {
    ($name:ident $(,)? $($init:ident),*) => {
        projection_delegate!{ $name, $name, $($init,)* }
    };
}

use downcast;
use projection;

const NUM_PROJECTIONS: usize = 17;

macro_rules! declare_projections {
    ($(($name:ident $(,)? $($init:ident),*)),+ $(,)?) => {
        const PROJECTIONS: [ProjInit; NUM_PROJECTIONS] = [
        $(
            ProjInit(stringify!($name), $name::stub::$name),
            $(
                ProjInit(stringify!($init), $name::stub::$init),
            )*
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

pub mod aea;
pub mod estmerc;
pub mod etmerc;
pub mod geocent;
pub mod laea;
pub mod latlong;
pub mod lcc;
pub mod merc;
pub mod somerc;
pub mod stere;
pub mod sterea;
pub mod tmerc;

#[rustfmt::skip]
declare_projections! [
    (latlong, longlat),
    (lcc),
    (etmerc, utm),
    (tmerc),
    (aea, leac),
    (stere, ups),
    (sterea),
    (merc, webmerc),
    (geocent, cart),
    (somerc),
    (laea),
];

///
/// Return the projection definition
///
pub(crate) fn find_projection(name: &str) -> Option<&ProjInit> {
    PROJECTIONS
        .iter()
        .find(|d| d.name().eq_ignore_ascii_case(name))
}
