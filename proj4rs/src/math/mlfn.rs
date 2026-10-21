//!
//! mlfn
//!  Meridional distance
//!
//! Meridional distance for ellipsoid and inverse using 6th-order expansion in
//! the third flattening n. This gives full double precision accuracy for |f| <= 1/150.
//!
//! Port of PROJ `src/mlfn.cpp` (PROJ >= 9.2, <https://github.com/OSGeo/PROJ/pull/3516>).
//!
use super::auxlat::{
    AuxCoeffs, C_MU_PHI, C_PHI_MU, coeffs_even, convert, rectifying_radius, third_flattening,
};

/// Meridional distance coefficients
#[derive(Debug, Clone, Copy)]
pub(crate) struct Enfn {
    /// Rectifying radius
    rr: f64,
    /// phi -> mu
    mu_phi: AuxCoeffs,
    /// mu -> phi
    phi_mu: AuxCoeffs,
}

/// Compute the meridional distance coefficients from the
/// first eccentricity squared.
pub(crate) fn enfn(es: f64) -> Enfn {
    let n = third_flattening(es);
    Enfn {
        rr: rectifying_radius(n),
        mu_phi: coeffs_even(n, &C_MU_PHI),
        phi_mu: coeffs_even(n, &C_PHI_MU),
    }
}

/// Meridional distance (in units of the semi-major axis)
#[inline]
pub(crate) fn mlfn(phi: f64, sphi: f64, cphi: f64, en: Enfn) -> f64 {
    en.rr * convert(phi, sphi, cphi, &en.mu_phi)
}

/// Inverse meridional distance: latitude from meridional distance
#[inline]
pub(crate) fn inv_mlfn(mu: f64, en: Enfn) -> f64 {
    let mu = mu / en.rr;
    convert(mu, mu.sin(), mu.cos(), &en.phi_mu)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::consts::FRAC_PI_2;

    #[test]
    fn mlfn_sphere() {
        let en = enfn(0.);
        for phi in [-1.5, -0.5, 0., 0.3, 1.2] {
            assert_eq!(mlfn(phi, phi.sin(), phi.cos(), en), phi);
            assert_eq!(inv_mlfn(phi, en), phi);
        }
    }

    #[test]
    fn mlfn_wgs84() {
        // WGS84
        let a = 6378137.;
        let rf = 298.257223563;
        let f = 1. / rf;
        let es = f * (2. - f);
        let en = enfn(es);

        // Quarter meridian for WGS84 (Karney, GeographicLib)
        let qm = a * mlfn(FRAC_PI_2, 1., 0., en);
        assert!((qm - 10001965.729312723).abs() < 1e-8, "{qm}");

        // Round trip
        for i in -90..=90 {
            let phi = (i as f64).to_radians();
            let mu = mlfn(phi, phi.sin(), phi.cos(), en);
            assert!((inv_mlfn(mu, en) - phi).abs() < 1e-15, "{i}");
        }
    }
}
