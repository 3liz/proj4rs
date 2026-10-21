use super::consts::EPS_7;

/// Computes coefficient q such that authalic latitude = beta = asin(q / qp)
/// where qp is q at phi=90deg, i.e. qp = qsfn(1, e, one_es)
///
/// Cf  Snyder (3-11) and (3-12)
///
/// From PROJ `pj_authalic_lat_q` (PROJ >= 9.7)
pub(crate) fn qsfn(sinphi: f64, e: f64, one_es: f64) -> f64 {
    if e >= EPS_7 {
        let e_sinphi = e * sinphi;
        let one_minus_e_sinphi_sq = 1.0 - e_sinphi * e_sinphi;
        // avoid zero division, fail gracefully
        if one_minus_e_sinphi_sq == 0.0 {
            f64::INFINITY
        } else {
            // Snyder uses 0.5 * ln((1-e*sinphi)/(1+e*sinphi) which is
            // -atanh(e*sinphi)
            one_es * (sinphi / one_minus_e_sinphi_sq + e_sinphi.atanh() / e)
        }
    } else {
        sinphi + sinphi
    }
}
