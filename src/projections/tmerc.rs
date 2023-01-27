//!
//! Transverse mercator
//!
//! Provide both Evenden/Snyder and Poder/ensager algorithm, depending on
//! parameters
//!
//! The default algorithm is Poder/Ensager except for the spherical case
//! where the Evenden/Snyder is used
//!
//!

use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::proj::ProjData;
use crate::projections::{estmerc, etmerc};

// Projection stub
super::projection! { tmerc }

#[derive(Debug)]
pub(crate) enum Projection {
    Exact(etmerc::Projection),
    Approx(estmerc::Projection),
}

use Projection::*;

impl Projection {
    const ALG_PARAM: &str = "algo";

    pub fn tmerc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        if p.ellps.is_sphere() || params.check_option("approx")? {
            Ok(Approx(estmerc::Projection::estmerc(p, params)?))
        } else {
            // try 'algo' parameter
            match params.try_value(Self::ALG_PARAM)? {
                Some("evenden_snyder") => Ok(Approx(estmerc::Projection::estmerc(p, params)?)),
                Some("poder_engsager") | None => Ok(Exact(etmerc::Projection::etmerc(p, params)?)),
                Some(_) => Err(Error::InvalidParameterValue(Self::ALG_PARAM)),
            }
        }
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Exact(p) => p.forward(lam, phi, z),
            Approx(p) => p.forward(lam, phi, z),
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Exact(p) => p.inverse(x, y, z),
            Approx(p) => p.inverse(x, y, z),
        }
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
    use crate::adaptors::transform_xy;
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};
    use approx::assert_abs_diff_eq;

    #[test]
    fn proj_estmerc_ell() {
        let p = Proj::from_proj_string("+proj=tmerc +ellps=GRS80 +approx").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222650.79679577847, 110642.22941192707, 0.)),
            ((2., -1., 0.), (222650.79679577847, -110642.22941192707, 0.)),
            ((-2., 1., 0.), (-222650.79679577847, 110642.22941192707, 0.)),
            (
                (-2., -1., 0.),
                (-222650.79679577847, -110642.22941192707, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_estmerc_sph() {
        // Spherical planet will choose estmerc algorithm
        let p = Proj::from_proj_string("+proj=tmerc +R=6400000").unwrap();

        println!("{:#?}", p.data());
        println!("{:#?}", p.projection());

        // Sames results as Proj9 'proj -d 11 +proj=tmerc +R=6400000  +approx'

        let inputs = [
            ((2., 1., 0.), (223413.46640632232, 111769.14504059685, 0.)),
            ((2., -1., 0.), (223413.46640632232, -111769.14504059685, 0.)),
            ((-2., 1., 0.), (-223413.46640632208, 111769.14504059685, 0.)),
            (
                (-2., -1., 0.),
                (-223413.46640632208, -111769.14504059685, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
