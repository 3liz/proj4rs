//!
//! Proj4 datum definitions
//!
use crate::ellipsoids::{constants::*, EllipsoidDefn};

/// Shift method is either
/// defined by Helmert transforms or nadgrids
pub enum ShiftMethod {
    ToWGS84_0,
    ToWGS84_3(f64, f64, f64),
    ToWGS84_7(f64, f64, f64, f64, f64, f64, f64),
    NadGrids(&'static str),
}

use ShiftMethod::*;

/// Static datums table
const DATUMS: [(&str, ShiftMethod, &EllipsoidDefn, &str); 17] = [
    ("WGS84", ToWGS84_0, &WGS84, ""),
    (
        "GGRS87",
        ToWGS84_3(-199.87, 74.79, 246.62),
        &GRS80,
        "Greek_Geodetic_Reference_System_1987",
    ),
    ("NAD83", ToWGS84_0, &GRS80, "North_American_Datum_1983"),
    (
        "NAD27",
        NadGrids("@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"),
        &CLRK66,
        "North_American_Datum_1927",
    ),
    // defn is "nadgrids=@BETA2007.gsb" in proj 9
    (
        "potsdam",
        ToWGS84_7(598.1, 73.7, 418.2, 0.202, 0.045, -2.455, 6.7),
        &BESSEL,
        "Potsdam Rauenberg 1950 DHDN",
    ),
    (
        "carthage",
        ToWGS84_3(-263.0, 6.0, 431.0),
        &CLRK80IGN,
        "Carthage 1934 Tunisia",
    ),
    (
        "hermannskogel",
        ToWGS84_7(577.326, 90.129, 463.919, 5.137, 1.474, 5.297, 2.4232),
        &BESSEL,
        "Hermannskogel",
    ),
    (
        "ire65",
        ToWGS84_7(482.530, -130.596, 564.557, -1.042, -0.214, -0.631, 8.15),
        &MOD_AIRY,
        "Ireland 1965",
    ),
    (
        "nzgd49",
        ToWGS84_7(59.47, -5.04, 187.44, 0.47, -0.1, 1.024, -4.5993),
        &INTL,
        "New Zealand Geodetic Datum 1949",
    ),
    (
        "OSGB36",
        ToWGS84_7(446.448, -125.157, 542.060, 0.1502, 0.2470, 0.8421, -20.4894),
        &AIRY,
        "Airy 1830",
    ),
    // Added from proj4js
    (
        "ch1903",
        ToWGS84_3(674.374, 15.056, 405.346),
        &BESSEL,
        "swiss",
    ),
    (
        "osni52",
        ToWGS84_7(482.530, -130.596, 564.557, -1.042, -0.214, -0.631, 8.15),
        &AIRY,
        "Irish National",
    ),
    (
        "rassadiran",
        ToWGS84_3(-133.63, -157.5, -158.62),
        &INTL,
        "Rassadiran",
    ),
    (
        "s_jtsk",
        ToWGS84_3(589., 76., 480.),
        &BESSEL,
        "S-JTSK (Ferro)",
    ),
    (
        "beduaram",
        ToWGS84_3(-106., -87., 188.),
        &CLRK80,
        "Beduaram",
    ),
    (
        "gunung_segara",
        ToWGS84_3(-403., 684., 41.),
        &BESSEL,
        "Gunung Segara Jakarta",
    ),
    (
        "rnb72",
        ToWGS84_7(106.869, -52.2978, 103.724, -0.33657, 0.456955, -1.84218, 1.),
        &INTL,
        "Reseau National Belge 1972",
    ),
];

/// Return the datum definition
pub fn datum_defn(name: &str) -> Option<(&ShiftMethod, &EllipsoidDefn)> {
    DATUMS
        .iter()
        .find(|d| d.0.eq_ignore_ascii_case(name))
        .map(|d| (&d.1, d.2))
}
