//!
//! Conversions between geographic and authalic latitude
//!
//! From PROJ `src/latitudes.cpp` (PROJ >= 9.7)
//!
use super::auxlat::{AuxCoeffs, C_PHI_XI, C_XI_PHI, coeffs_full, convert, third_flattening};
use super::qsfn;
use crate::ellps::Ellipsoid;

/// Authalic latitude conversions for a given ellipsoid
#[derive(Debug, Clone)]
pub(crate) struct Authalic {
    e: f64,
    es: f64,
    one_es: f64,
    /// q at phi=90deg
    qp: f64,
    /// Use the series in n for the forward conversion and
    /// skip the Newton refinement of the inverse conversion
    series: bool,
    /// xi -> phi
    phi_xi: AuxCoeffs,
    /// phi -> xi
    xi_phi: AuxCoeffs,
}

impl Authalic {
    pub fn new(ellps: &Ellipsoid) -> Self {
        let (e, es, one_es) = (ellps.e, ellps.es, ellps.one_es);
        let n = third_flattening(es);
        // The series are preferred over the analytical expressions for
        // |n| < 0.01.  However the inverse series is used to start the
        // inverse method for large |n|.
        let series = n.abs() < 0.01;
        Self {
            e,
            es,
            one_es,
            qp: qsfn(1., e, one_es),
            series,
            phi_xi: coeffs_full(n, &C_PHI_XI),
            xi_phi: if series {
                coeffs_full(n, &C_XI_PHI)
            } else {
                [0.; _]
            },
        }
    }

    /// q at phi=90deg
    #[inline]
    pub fn qp(&self) -> f64 {
        self.qp
    }

    /// Coefficient q such that authalic latitude = asin(q / qp)
    #[inline]
    pub fn q(&self, sinphi: f64) -> f64 {
        qsfn(sinphi, self.e, self.one_es)
    }

    /// Authalic latitude from the geographic latitude
    pub fn lat(&self, phi: f64, sinphi: f64, cosphi: f64) -> f64 {
        if self.series {
            convert(phi, sinphi, cosphi, &self.xi_phi)
        } else {
            // This result is ill-conditioned near the poles.  So don't use this if
            // the series are accurate.
            (self.q(sinphi) / self.qp).clamp(-1., 1.).asin()
        }
    }

    /// Geographic latitude from the authalic latitude
    pub fn lat_inverse(&self, beta: f64) -> f64 {
        let mut phi = convert(beta, beta.sin(), beta.cos(), &self.phi_xi);
        if self.series {
            return phi;
        }
        // If the flattening is large, solve
        //   f(phi) = qp*sin(beta)/(1-e^2) - q(phi)/(1-e^2) = 0
        // for phi, using Newton's method, where
        //   q(phi)/(1-e^2) = sin(phi)/(1 - e^2*sin(phi)^2) + atanh(e*sin(phi))/e
        // and
        //   df(phi)/dphi = - dq(phi)/dphi / (1-e^2)
        //                = - 2 * (1-e^2) * cos(phi) / (1 - e^2*sinphi^2)^2
        // This is subject to large roundoff errors near the poles, so only use
        // this if the series isn't accurate.
        let q_div_one_minus_es = beta.sin() * self.qp / self.one_es;
        for _ in 0..10 {
            let (sinphi, cosphi) = phi.sin_cos();
            let one_minus_es_sin2phi = 1. - self.es * (sinphi * sinphi);
            let dphi = (one_minus_es_sin2phi * one_minus_es_sin2phi) / (2. * cosphi)
                * (q_div_one_minus_es
                    - sinphi / one_minus_es_sin2phi
                    - (self.e * sinphi).atanh() / self.e);
            if dphi.is_nan() || dphi.abs() < 1e-15 {
                break;
            }
            phi += dphi;
        }
        phi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trip(ellps: &Ellipsoid, max_lat: i32, prec: f64) {
        let auth = Authalic::new(ellps);
        for i in -max_lat..=max_lat {
            let phi = (i as f64).to_radians();
            let xi = auth.lat(phi, phi.sin(), phi.cos());
            assert!((auth.lat_inverse(xi) - phi).abs() < prec, "{i}");
        }
    }

    #[test]
    fn authalic_round_trip_grs80() {
        let ellps = Ellipsoid::try_from_ellipsoid(&crate::ellipsoids::constants::GRS80).unwrap();
        assert!(Authalic::new(&ellps).series);
        round_trip(&ellps, 90, 1e-15);
    }

    #[test]
    fn authalic_round_trip_large_flattening() {
        // |n| > 0.01: analytical forward and Newton inverse
        let ellps = Ellipsoid::calc_ellipsoid_params(1., crate::ellps::Shape::SP_rf(10.)).unwrap();
        assert!(!Authalic::new(&ellps).series);
        // As in PROJ, the Newton refinement is ill-conditioned at the poles
        round_trip(&ellps, 89, 1e-12);
    }
}
