//!
//! Determine latitude from authalic latitude
//!

pub(crate) fn authset(es: f64) -> (f64, f64, f64) {
    const P00: f64 = 1. / 3.;
    const P01: f64 = 31. / 180.;
    const P02: f64 = 517. / 5040.;
    const P10: f64 = 23. / 360.;
    const P11: f64 = 251. / 3780.;
    const P20: f64 = 761. / 45360.;
    let t = es * es;
    (
        es * P00 + t * P01 + t * es * P02,
        t * P10 + t * es * P11,
        t * es * P20,
    )
}

pub(crate) fn authlat(beta: f64, apa: (f64, f64, f64)) -> f64 {
    let t = beta + beta;
    beta + apa.0 * t.sin() + apa.1 * (t + t).sin() + apa.2 * (t + t + t).sin()
}
