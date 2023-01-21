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

// Projection stub
super::projection!(latlong);

pub(super) const NAME: &str = "latlon";

#[derive(Debug)]
pub(crate) struct Projection {}

impl Projection {
    pub fn init(p: &mut Proj, params: &ParamList) -> Result<Self> {
        p.is_latlong = true;
        p.x0 = Some(0.);
        p.y0 = Some(0.);
        Ok(Self {})
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam, phi, z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((x, y, z))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adaptors::transform_xy;
    use crate::proj::Proj;

    #[test]
    fn proj_latlon_init() {
        let p = Proj::from_proj_string("+proj=latlon +datum=WGS84").unwrap();

        assert_eq!(p.x0.unwrap(), 0.);
        assert_eq!(p.y0.unwrap(), 0.);
        assert_eq!(p.projname(), NAME);
    }

    #[test]
    fn proj_latlon_to_latlon() {
        let p_from = Proj::from_proj_string("+proj=latlon +datum=WGS84").unwrap();
        let p_to = Proj::from_proj_string("+proj=latlon +datum=WGS84").unwrap();

        let (lon_in, lat_in) = (2.3522219, 48.856614);

        let (lon_out, lat_out) = transform_xy(&p_from, &p_to, lon_in, lat_in).unwrap();
        assert_eq!((lon_out, lat_out), (lon_in, lat_in));
    }
}
