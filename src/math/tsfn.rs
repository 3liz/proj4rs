use super::consts::FRAC_PI_2;

#[inline]
pub(crate) fn tsfn(phi: f64, sinphi: f64, e: f64) -> f64 {
    //  XXX Avoid division by zero, check denominator
    (0.5 * (FRAC_PI_2 - phi)).tan() / ((1. - sinphi * e) / (1. + sinphi * e)).powf(0.5 * e)
}
