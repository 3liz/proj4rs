//!
//! From proj/eqc.cpp
//!
//! See also <https://proj.org/operations/projections/eqc.html>
//!
//! Simplest of all projections
//!
//! eqc: "Equidistant Cylindrical (Plate Carree)"
//!

use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { eqc }

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    rc: f64,
    phi0: f64,
}

impl Projection {
    pub fn eqc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let rc = params.try_angular_value("lat_ts")?.unwrap_or(0.).cos();
        if rc <= 0. {
            return Err(Error::InvalidParameterValue("lat_ts should be <= 90°"));
        }
        p.ellps = Ellipsoid::sphere(p.ellps.a)?;
        Ok(Self { rc, phi0: p.phi0 })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((lam * self.rc, phi - self.phi0, z))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok((x / self.rc, y + self.phi0, z))
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
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_eqc_wgs84() {
        let p = Proj::from_proj_string("+proj=eqc +ellps=WGS84").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [((2., 47., 0.), (222638.98158654713, 5232016.06728385761, 0.))];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_eqc_lat_ts() {
        let p = Proj::from_proj_string("+proj=eqc +lat_ts=30 +lon_0=-90").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (-88., 30., 0.),
            (192811.01392664597, 3339584.72379820701, 0.),
        )];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
