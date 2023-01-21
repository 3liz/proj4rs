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

pub(crate) type ProjFn = fn(&Proj, f64, f64, f64) -> Result<(f64, f64, f64)>;

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

macro_rules! proj {
    ($name:ident) => {
        ProjInit($name::NAME, $name::Projection::init)
    };
}

// Macro for retrieval of parameters from the proj object
// not that makes us writing a match to a unique element.
macro_rules! downcast {
    ($name:ident, $p:expr) => {
        match &$p.projdata {
            crate::projections::ProjParams::$name(data) => data,
            _ => unreachable!(),
        }
    };
}

pub(crate) use downcast;

/// Hold per projection data calculated in the `init`
/// function
#[allow(non_camel_case_types)]
#[derive(Debug)]
pub(crate) enum ProjParams {
    NoParams,
    latlong(latlong::Projection),
    lcc(lcc::Projection),
}

// ----------------------------
// Projection list
// ---------------------------
const PROJECTIONS: [ProjInit; 2] = [proj!(latlong), proj!(lcc)];

/// Return the datum definition
pub(crate) fn find_projection(name: &str) -> Option<&ProjInit> {
    PROJECTIONS
        .iter()
        .find(|d| d.name().eq_ignore_ascii_case(name))
}
