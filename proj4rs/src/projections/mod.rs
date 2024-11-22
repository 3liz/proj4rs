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
#[derive(Clone)]
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
    #[inline(always)]
    pub fn has_inverse(&self) -> bool {
        self.3
    }
    #[inline(always)]
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

macro_rules! declare_proj {
    ($name:ident) => {
        ProjInit(stringify!($name), $name::stub::$name)
    };
    ($name:ident, $init:ident) => {
        ProjInit(stringify!($init), $name::stub::$init)
    };
}

// ----------------------------
// Projection list
// ---------------------------

pub mod aea;
pub mod eqc;
pub mod estmerc;
pub mod etmerc;
pub mod geocent;
pub mod geos;
pub mod laea;
pub mod latlong;
pub mod lcc;
pub mod merc;
pub mod moll;
pub mod somerc;
pub mod stere;
pub mod sterea;
pub mod tmerc;

#[cfg(feature = "aeqd")]
pub mod aeqd;

#[allow(unused_mut)]
const fn num_projections() -> usize {
    let mut num = 22;
    #[cfg(feature = "aeqd")]
    {
        num += 1;
    }
    num
}

const NUM_PROJECTIONS: usize = num_projections();

#[rustfmt::skip]
const PROJECTIONS: [ProjInit; NUM_PROJECTIONS] = [
    declare_proj!(latlong),
    declare_proj!(latlong, longlat),
    declare_proj!(lcc),
    declare_proj!(etmerc),
    declare_proj!(etmerc, utm),
    declare_proj!(tmerc),
    declare_proj!(aea),
    declare_proj!(aea, leac),
    declare_proj!(stere),
    declare_proj!(stere, ups),
    declare_proj!(sterea),
    declare_proj!(merc),
    declare_proj!(merc, webmerc),
    declare_proj!(geocent),
    declare_proj!(geocent, cart),
    declare_proj!(somerc),
    declare_proj!(laea),
    declare_proj!(moll),
    declare_proj!(moll, wag4),
    declare_proj!(moll, wag5),
    declare_proj!(geos),
    declare_proj!(eqc),
    #[cfg(feature = "aeqd")]
    declare_proj!(aeqd),
];

#[allow(non_camel_case_types)]
#[derive(Debug, Clone)]
pub(crate) enum ProjParams {
    latlong(latlong::Projection),
    lcc(lcc::Projection),
    etmerc(etmerc::Projection),
    tmerc(tmerc::Projection),
    aea(aea::Projection),
    stere(stere::Projection),
    sterea(sterea::Projection),
    merc(merc::Projection),
    geocent(geocent::Projection),
    somerc(somerc::Projection),
    laea(laea::Projection),
    moll(moll::Projection),
    geos(geos::Projection),
    eqc(eqc::Projection),
    #[cfg(feature = "aeqd")]
    aeqd(aeqd::Projection),
}

///
/// Return the projection definition
///
pub(crate) fn find_projection(name: &str) -> Option<&ProjInit> {
    PROJECTIONS
        .iter()
        .find(|d| d.name().eq_ignore_ascii_case(name))
}
