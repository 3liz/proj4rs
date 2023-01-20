//!
//! Most important projections:
//!
//! lcc, merc, tmerc, utm et aea
//! stere et sterea pour les poles.
//!

mod latlong;

use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::Proj;

pub(crate) type InitFn = fn(&mut Proj, &ParamList) -> Result<ProjParams>;

pub(crate) type ProjFn = fn(&Proj, f64, f64, f64) -> Result<(f64, f64, f64)>;

struct ProjDecl(&'static str, InitFn, ProjFn, ProjFn);

/// Hold per projection data calculated in the `init`
/// function
#[allow(non_camel_case_types)]
#[derive(Debug)]
pub(crate) enum ProjParams {
    NoParams,
    latlong(latlong::Projection),
}

macro_rules! proj {
    ($name:ident) => {
        ProjDecl(
            $name::NAME,
            $name::Projection::init,
            $name::Projection::inverse,
            $name::Projection::forward,
        )
    };
}

const PROJECTIONS: [ProjDecl; 1] = [proj!(latlong)];
