//!
//! From proj/eqc.cpp
//!
//! See also <https://proj.org/operations/projections/eqc.html>
//!
//! Spherical (EPSG method 1029) and, with the default `epsg-method-1028` feature,
//! ellipsoidal (EPSG method 1028) forms, as in PROJ >= 9.8
//! (<https://github.com/OSGeo/PROJ/pull/4656>)
//!
//! eqc: "Equidistant Cylindrical (Plate Carree)"
//!

#[cfg(not(feature = "epsg-method-1028"))]
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
#[cfg(feature = "epsg-method-1028")]
use crate::math::{Enfn, enfn, inv_mlfn, mlfn};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { eqc }

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    /// Spherical (EPSG method 1029)
    Sph { rc: f64, phi0: f64 },
    /// Ellipsoidal (EPSG method 1028)
    #[cfg(feature = "epsg-method-1028")]
    Ell { rc: f64, m0: f64, en: Enfn },
}

impl Projection {
    pub fn eqc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let phi1 = params.try_angular_value("lat_ts")?.unwrap_or(0.);
        if phi1.cos() <= 0. {
            return Err(Error::InvalidParameterValue("lat_ts should be <= 90°"));
        }

        #[cfg(feature = "epsg-method-1028")]
        if p.ellps.is_ellipsoid() {
            return Ok(Self::ellipsoidal(p, phi1));
        }

        // Without EPSG method 1028, use the spherical formula
        // on a sphere of radius a, as PROJ < 9.8
        #[cfg(not(feature = "epsg-method-1028"))]
        {
            p.ellps = Ellipsoid::sphere(p.ellps.a)?;
        }

        Ok(Self::Sph {
            rc: phi1.cos(),
            phi0: p.phi0,
        })
    }

    #[cfg(feature = "epsg-method-1028")]
    fn ellipsoidal(p: &ProjData, phi1: f64) -> Self {
        let es = p.ellps.es;
        let (sin_phi1, cos_phi1) = phi1.sin_cos();
        // nu1 = a / sqrt(1 - e^2 sin^2(lat_ts)), normalized by a
        let nu1 = 1. / (1. - es * sin_phi1 * sin_phi1).sqrt();
        let en = enfn(es);
        Self::Ell {
            rc: nu1 * cos_phi1,
            m0: mlfn(p.phi0, p.phi0.sin(), p.phi0.cos(), en),
            en,
        }
    }

    // Ellipsoidal (EPSG method 1028):
    //   E = FE + nu1 cos(lat_ts) (lon - lon_0)
    //   N = FN + M - M0
    // where M is the meridional arc from the equator to the latitude
    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok(match self {
            Self::Sph { rc, phi0 } => (lam * rc, phi - phi0, z),
            #[cfg(feature = "epsg-method-1028")]
            Self::Ell { rc, m0, en } => (lam * rc, mlfn(phi, phi.sin(), phi.cos(), *en) - m0, z),
        })
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        Ok(match self {
            Self::Sph { rc, phi0 } => (x / rc, y + phi0, z),
            #[cfg(feature = "epsg-method-1028")]
            Self::Ell { rc, m0, en } => (x / rc, inv_mlfn(y + m0, *en), z),
        })
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

    // Ellipsoidal reference values: x from the EPSG 1028 formula,
    // y = meridian arc, from PROJ 9.4 `+proj=tmerc` on its central meridian

    #[test]
    #[cfg(feature = "epsg-method-1028")]
    fn proj_eqc_wgs84() {
        let p = Proj::from_proj_string("+proj=eqc +ellps=WGS84").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [((2., 47., 0.), (222638.981586547, 5207247.008955783, 0.))];

        test_proj_forward(&p, &inputs, 1e-9);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    #[cfg(feature = "epsg-method-1028")]
    fn proj_eqc_lat_ts() {
        // Default ellipsoid (WGS84)
        let p = Proj::from_proj_string("+proj=eqc +lat_ts=30 +lon_0=-90").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [((-88., 30., 0.), (192972.560501793, 3320113.397940382, 0.))];

        test_proj_forward(&p, &inputs, 1e-9);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    #[cfg(feature = "epsg-method-1028")]
    fn proj_eqc_lat_0() {
        let p = Proj::from_proj_string("+proj=eqc +ellps=WGS84 +lat_0=45").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((0., 45., 0.), (0., 0., 0.)),
            ((0., 60., 0.), (0., 1669128.441512768, 0.)),
            ((2., 47., 0.), (222638.981586547, 222302.630978040, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-9);
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_eqc_sphere() {
        // Spherical formulation (EPSG method 1029), values from PROJ 9.4
        let p = Proj::from_proj_string("+proj=eqc +lat_ts=30 +lon_0=-90 +R=6378137").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (-88., 30., 0.),
            (192811.01392664597, 3339584.72379820701, 0.),
        )];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    // Without EPSG method 1028: spherical formula on a sphere of radius a
    // values from PROJ 9.4

    #[test]
    #[cfg(not(feature = "epsg-method-1028"))]
    fn proj_eqc_wgs84_spherical() {
        let p = Proj::from_proj_string("+proj=eqc +ellps=WGS84").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [((2., 47., 0.), (222638.98158654713, 5232016.06728385761, 0.))];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    #[cfg(not(feature = "epsg-method-1028"))]
    fn proj_eqc_lat_ts_spherical() {
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
