//!
//! Stub projection implementation for lat/long coordinates.
//!
//! We don't actually change the coordinates, but we want proj=latlong
//! to act sort of like a projection.
//!
use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::{ProjData, ProjType};

// Projection stub
super::projection! { latlong, longlat }

#[derive(Debug)]
pub(crate) struct Projection {}

impl Projection {
    pub fn latlong(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        p.proj_type = ProjType::Latlong;
        p.x0 = 0.;
        p.y0 = 0.;
        Ok(Self {})
    }

    pub fn longlat(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::latlong(p, params)
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam, phi, z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((x, y, z))
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::adaptors::transform_xy;
    use crate::proj::Proj;

    #[test]
    fn proj_latlon_init() {
        let p = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

        let d = p.data();

        assert_eq!(d.x0, 0.);
        assert_eq!(d.y0, 0.);
        assert_eq!(p.projname(), "latlong");
    }

    #[test]
    fn proj_latlon_to_latlon() {
        let p_from = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();
        let p_to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

        let (lon_in, lat_in) = (2.3522219, 48.856614);

        let (lon_out, lat_out) = transform_xy(&p_from, &p_to, lon_in, lat_in).unwrap();
        assert_eq!((lon_out, lat_out), (lon_in, lat_in));
    }
}
