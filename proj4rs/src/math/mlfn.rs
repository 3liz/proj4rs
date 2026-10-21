//!
//! mlfn
//!  Meridional distance
//!
//! Meridional distance for ellipsoid and inverse using 6th-order expansion in
//! the third flattening n. This gives full double precision accuracy for |f| <= 1/150.
//!
//! Port of PROJ `src/mlfn.cpp` (PROJ >= 9.2, <https://github.com/OSGeo/PROJ/pull/3516>).
//! See C. F. F. Karney, "On auxiliary latitudes", <https://arxiv.org/abs/2212.05818>
//!

const LMAX: usize = 6;

/// Meridional distance coefficients:
/// [0]: overall multiplier, [1..=LMAX]: phi -> mu, [LMAX+1..=2*LMAX]: mu -> phi
pub(crate) type Enfn = [f64; 2 * LMAX + 1];

/// Evaluate sum(p[i] * x^i, i, 0, N) via Horner's method.
#[inline]
fn polyval(x: f64, p: &[f64]) -> f64 {
    p.iter().rev().fold(0., |y, c| y * x + c)
}

/// Evaluate y = sum(c[k] * sin((2*k+2) * zeta), k, 0, K-1)
#[inline]
fn clenshaw(szeta: f64, czeta: f64, c: &[f64]) -> f64 {
    // Approx operation count = (K + 5) mult and (2 * K + 2) add
    let x = 2. * (czeta - szeta) * (czeta + szeta); // 2 * cos(2*zeta)
    let (u0, _) = c.iter().rev().fold((0., 0.), |(u0, u1), ck| (x * u0 - u1 + ck, u0));
    2. * szeta * czeta * u0 // sin(2*zeta) * u0
}

/// Compute the meridional distance coefficients from the
/// first eccentricity squared.
pub(crate) fn enfn(es: f64) -> Enfn {
    // Expansion of (quarter meridian) / ((a+b)/2 * pi/2) as series in n^2;
    // these coefficients are ( (2*k - 3)!! / (2*k)!! )^2 for k = 0..3
    const COEFF_RAD: [f64; 4] = [1., 1. / 4., 1. / 64., 1. / 256.];

    // Coefficients to convert phi to mu, Eq. A5 in arXiv:2212.05818
    // with 0 terms dropped
    #[rustfmt::skip]
    const COEFF_MU_PHI: [f64; 12] = [
        -3. / 2.,    9. / 16.,     -3. / 32.,     15. / 16.,
        -15. / 32.,  135. / 2048., -35. / 48.,    105. / 256.,
        315. / 512., -189. / 512., -693. / 1280., 1001. / 2048.,
    ];
    // Coefficients to convert mu to phi, Eq. A6 in arXiv:2212.05818
    // with 0 terms dropped
    #[rustfmt::skip]
    const COEFF_PHI_MU: [f64; 12] = [
        3. / 2.,      -27. / 32.,      269. / 512.,   21. / 16.,
        -55. / 32.,   6759. / 4096.,   151. / 96.,    -417. / 128.,
        1097. / 512., -15543. / 2560., 8011. / 2560., 293393. / 61440.,
    ];

    // Third flattening: n = (1 - b/a) / (1 + b/a) = es / (1 + sqrt(1 - es))^2
    let n = es / (1. + (1. - es).sqrt()).powi(2);
    let n2 = n * n;

    let mut en = [0.; 2 * LMAX + 1];
    en[0] = polyval(n2, &COEFF_RAD[..=LMAX / 2]) / (1. + n);

    let mut d = n;
    let mut o = 0;
    for l in 0..LMAX {
        let m = (LMAX - l - 1) / 2;
        en[l + 1] = d * polyval(n2, &COEFF_MU_PHI[o..=o + m]);
        en[l + 1 + LMAX] = d * polyval(n2, &COEFF_PHI_MU[o..=o + m]);
        d *= n;
        o += m + 1;
    }
    en
}

/// Meridional distance (in units of the semi-major axis)
#[inline]
pub(crate) fn mlfn(phi: f64, sphi: f64, cphi: f64, en: Enfn) -> f64 {
    en[0] * (phi + clenshaw(sphi, cphi, &en[1..=LMAX]))
}

/// Inverse meridional distance: latitude from meridional distance
#[inline]
pub(crate) fn inv_mlfn(mu: f64, en: Enfn) -> f64 {
    let mu = mu / en[0];
    mu + clenshaw(mu.sin(), mu.cos(), &en[LMAX + 1..])
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
