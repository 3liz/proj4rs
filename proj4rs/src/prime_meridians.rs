///
/// Static prime meridians definitions
///
#[rustfmt::skip]
pub(crate) const PRIME_MERIDIANS: [(&str, &str, f64); 14] = [
    ("greenwich",  "0dE",              0.0),
    ("lisbon",     "9d07'54.862\"W",   -9.131906111111),
    ("paris",      "2d20'14.025\"E",   2.337229166667),
    ("bogota",     "74d04'51.3\"W",    -74.080916666667),
    ("madrid",     "3d41'16.58\"W",    -3.687938888889),
    ("rome",       "12d27'8.4\"E",     12.452333333333),
    ("bern",       "7d26'22.5\"E",     7.439583333333),
    ("jakarta",    "106d48'27.79\"E",  106.807719444444),
    ("ferro",      "17d40'W",          -17.666666666667),
    ("brussels",   "4d22'4.71\"E",     4.367975),
    ("stockholm",  "18d3'29.8\"E",     18.058277777778),
    ("athens",     "23d42'58.815\"E",  23.7163375),
    ("oslo",       "10d43'22.5\"E",    10.722916666667),
    ("copenhagen", "12d34'40.35\"E",   12.57788),
];

/// Return the datum definition
pub fn find_prime_meridian(name: &str) -> Option<f64> {
    PRIME_MERIDIANS
        .iter()
        .find(|d| d.0.eq_ignore_ascii_case(name))
        .map(|d| d.2)
}
