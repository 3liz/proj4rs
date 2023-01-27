use super::consts::EPS_7;

pub(crate) fn qsfn(sinphi: f64, e: f64, one_es: f64) -> f64 {
    if e >= EPS_7 {
        let con = e * sinphi;
        let div1 = 1.0 - con * con;
        let div2 = 1.0 + con;
        // avoid zero division, fail gracefully
        if div1 == 0.0 || div2 == 0.0 {
            f64::INFINITY
        } else {
            one_es * (sinphi / div1 - (0.5 / e) * ((1. - con) / div2).ln())
        }
    } else {
        // XXX why not 2.*sinphi ?
        sinphi + sinphi
    }
}
