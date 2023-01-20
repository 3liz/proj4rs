//
// Stub projection implementation for lat/long coordinates. We
// don't actually change the coordinates, but we want proj=latlong
// to act sort of like a projection.
//
// Original author: Frank Warmerdam, warmerdam@pobox.com
//
use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::Proj;

use super::ProjParams;

pub(super) const NAME: &str = "latlong";

#[derive(Debug, Default)]
pub(crate) struct Projection {}

impl Projection {
    pub fn init(p: &mut Proj, params: &ParamList) -> Result<ProjParams> {
        p.is_latlong = true;
        p.x0 = 0.;
        p.y0 = 0.;
        p.inverse = Some(Self::inverse);
        p.forward = Some(Self::forward);
        Ok(ProjParams::latlong(Self {}))
    }

    pub fn forward(p: &Proj, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam / p.ellps.a, phi / p.ellps.a, z))
    }

    pub fn inverse(p: &Proj, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam * p.ellps.a, phi * p.ellps.a, z))
    }
}
