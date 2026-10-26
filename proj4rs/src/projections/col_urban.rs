//! Colombia Urban
//!
//! ref: <https://proj.org/operations/projections/col_urban.html>
//!
//! Miscellaneous, ellipsoidal, for Colombian cities. Takes `h_0`, the height of the
//! projection plane, which defaults to 0.
//!
//! Notations and formulas taken from IOGP Publication 373-7-2 -
//! Geomatics Guidance Note number 7, part 2 - March 2020
//!
use crate::errors::Result;
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { col_urban }

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    es: f64,
    phi0: f64,
    h0: f64,   // height of projection origin, divided by semi-major axis (a)
    rho0: f64, // adimensional value, contrary to Guidance note 7.2
    a: f64,
    b: f64, // adimensional value, contrary to Guidance note 7.2
    c: f64,
    d: f64, // adimensional value, contrary to Guidance note 7.2
}

impl Projection {
    pub fn col_urban(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let el = &p.ellps;
        let h0 = params.try_value::<f64>("h_0")?.unwrap_or(0.) / el.a;
        let sinphi0 = p.phi0.sin();
        let nu0 = 1. / (1. - el.es * sinphi0 * sinphi0).sqrt();
        let rho0 = (1. - el.es) / (1. - el.es * sinphi0 * sinphi0).powf(1.5);
        Ok(Self {
            es: el.es,
            phi0: p.phi0,
            h0,
            rho0,
            a: 1. + h0 / nu0,
            b: p.phi0.tan() / (2. * rho0 * nu0),
            c: 1. + h0,
            d: rho0 * (1. + h0 / (1. - el.es)),
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let nu = 1. / (1. - self.es * sinphi * sinphi).sqrt();
        let lam_nu_cosphi = lam * nu * cosphi;
        let sinphi_m = (0.5 * (phi + self.phi0)).sin();
        let rho_m = (1. - self.es) / (1. - self.es * sinphi_m * sinphi_m).powf(1.5);
        let g = 1. + self.h0 / rho_m;
        Ok((
            self.a * lam_nu_cosphi,
            g * self.rho0 * ((phi - self.phi0) + self.b * lam_nu_cosphi * lam_nu_cosphi),
            z,
        ))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let phi = self.phi0 + y / self.d - self.b * (x / self.c) * (x / self.c);
        let sinphi = phi.sin();
        let nu = 1. / (1. - self.es * sinphi * sinphi).sqrt();
        Ok((x / (self.c * nu * phi.cos()), phi, z))
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

//============
// Tests
//============

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Test values from proj builtins.gie
    #[test]
    fn proj_col_urban() {
        let p = Proj::from_proj_string(concat!(
            "+proj=col_urban +lat_0=4.68048611111111 +lon_0=-74.1465916666667",
            " +x_0=92334.879 +y_0=109320.965 +h_0=2550 +ellps=GRS80"
        ))
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [((-74.25, 4.8, 0.), (80859.033, 122543.174, 0.))];
        test_proj_forward(&p, &inputs, 1e-3);
        test_proj_inverse(&p, &inputs, 1e-7);
    }
}
