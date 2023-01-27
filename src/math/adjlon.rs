use super::consts::{EPS_12, PI, TAU};

pub(crate) fn adjlon(mut lon: f64) -> f64 {
    // Let lon slightly overshoot,
    // to avoid spurious sign switching at the date line
    if lon.abs() >= PI + EPS_12 {
        // adjust to 0..2pi rad
        lon += PI;

        // remove integral # of 'revolutions'
        lon -= TAU * (lon / TAU).floor();

        // adjust back to -pi..pi rad
        lon -= PI;
    }
    lon
}
