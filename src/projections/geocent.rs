//
// Stub projection implementation for geocent coordinates. No
// transformation occurs here because it is handled in transform.rs
//
use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::{ProjData, ProjType};

// Projection stub
super::projection! { geocent, cart }

#[derive(Debug)]
pub(crate) struct Projection {}

impl Projection {
    pub fn geocent(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        p.proj_type = ProjType::Geocentric;
        p.x0 = 0.;
        p.y0 = 0.;
        Ok(Self {})
    }

    pub fn cart(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        Self::geocent(p, params)
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
    use super::*;
    use crate::adaptors::transform_xyz;
    use crate::proj::Proj;
    use approx::assert_abs_diff_eq;

    #[test]
    fn proj_geocent() {
        let p_from = Proj::from_proj_string("+proj=latlong").unwrap();
        let p_to = Proj::from_proj_string("+proj=geocent +R=1").unwrap();

        let (lon_in, lat_in, z_in) = (0.0f64, 0.0f64, 0.0f64);

        let (x, y, z) = transform_xyz(
            &p_from,
            &p_to,
            lon_in.to_radians(),
            lat_in.to_radians(),
            z_in,
        )
        .unwrap();

        assert_abs_diff_eq!(x, 1.0, epsilon = 1.0e-8);
        assert_abs_diff_eq!(y, 0.0, epsilon = 1.0e-8);
        assert_abs_diff_eq!(z, 0.0, epsilon = 1.0e-8);
    }
}
