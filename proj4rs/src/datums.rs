//!
//! Proj4 datum definitions
//!
use crate::ellipsoids::{constants as ellps, EllipsoidDefn};

/// Shift method is either
/// defined by Helmert transforms or nadgrids
pub enum DatumParamDefn {
    ToWGS84_0,
    ToWGS84_3(f64, f64, f64),
    ToWGS84_7(f64, f64, f64, f64, f64, f64, f64),
    NadGrids(&'static str),
}

pub struct DatumDefn {
    pub id: &'static str,
    pub params: DatumParamDefn,
    pub ellps: &'static EllipsoidDefn,
    //pub comment: &'static str,
}

//#[rustfmt::skip]
pub mod constants {
    use super::*;

    macro_rules! nadgrids {
        ($grids:expr) => {
            DatumParamDefn::NadGrids($grids)
        };
    }

    macro_rules! towgs84 {
        ($x:expr, $y:expr, $z:expr) => {
            DatumParamDefn::ToWGS84_3($x, $y, $z)
        };
        ($x:expr, $y:expr, $z:expr, $rx:expr, $ry:expr, $rz:expr, $s:expr) => {
            DatumParamDefn::ToWGS84_7($x, $y, $z, $rx, $ry, $rz, $s)
        };
        () => {
            DatumParamDefn::ToWGS84_0
        };
    }

    macro_rules! datum {
        ($name:ident, $id:expr, $params:expr, $ellps:ident, $c:expr $(,)?) => {
            pub(crate) const $name: DatumDefn = DatumDefn {
                id: $id,
                params: $params,
                ellps: &ellps::$ellps,
                //comment: $c,
            };
        };
    }

    // ---------------------------
    //
    // Datum definitions
    //
    // ---------------------------
    datum!(WGS84, "WGS84", towgs84!(), WGS84, "");
    datum!(
        GGRS87,
        "GGRS87",
        towgs84!(-199.87, 74.79, 246.62),
        GRS80,
        "Greek_Geodetic_Reference_System_1987",
    );
    datum!(
        NAD83,
        "NAD83",
        towgs84!(),
        GRS80,
        "North_American_Datum_1983"
    );
    datum!(
        NAD27,
        "NAD27",
        nadgrids!("@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"),
        CLRK66,
        "North_American_Datum_1927",
    );
    // defn is "nadgrids=@BETA2007.gsb" in proj 9
    datum!(
        POTSDAM,
        "potsdam",
        towgs84!(598.1, 73.7, 418.2, 0.202, 0.045, -2.455, 6.7),
        BESSEL,
        "Potsdam Rauenberg 1950 DHDN",
    );
    datum!(
        CARTHAGE,
        "carthage",
        towgs84!(-263.0, 6.0, 431.0),
        CLRK80IGN,
        "Carthage 1934 Tunisia",
    );
    datum!(
        HERMANNSKOGEL,
        "hermannskogel",
        towgs84!(577.326, 90.129, 463.919, 5.137, 1.474, 5.297, 2.4232),
        BESSEL,
        "Hermannskogel",
    );
    datum!(
        IRE65,
        "ire65",
        towgs84!(482.530, -130.596, 564.557, -1.042, -0.214, -0.631, 8.15),
        MOD_AIRY,
        "Ireland 1965",
    );
    datum!(
        NZGD49,
        "nzgd49",
        towgs84!(59.47, -5.04, 187.44, 0.47, -0.1, 1.024, -4.5993),
        INTL,
        "New Zealand Geodetic Datum 1949",
    );
    datum!(
        OSGB36,
        "OSGB36",
        towgs84!(446.448, -125.157, 542.060, 0.1502, 0.2470, 0.8421, -20.4894),
        AIRY,
        "Airy 1830",
    );
    // Added from proj4js
    datum!(
        CH1903,
        "ch1903",
        towgs84!(674.374, 15.056, 405.346),
        BESSEL,
        "swiss",
    );
    datum!(
        OSNI52,
        "osni52",
        towgs84!(482.530, -130.596, 564.557, -1.042, -0.214, -0.631, 8.15),
        AIRY,
        "Irish National",
    );
    datum!(
        RASSADIRAN,
        "rassadiran",
        towgs84!(-133.63, -157.5, -158.62),
        INTL,
        "Rassadiran",
    );
    datum!(
        S_JTSK,
        "s_jtsk",
        towgs84!(589., 76., 480.),
        BESSEL,
        "S-JTSK (Ferro)",
    );
    datum!(
        BEDUARAM,
        "beduaram",
        towgs84!(-106., -87., 188.),
        CLRK80,
        "Beduaram",
    );
    datum!(
        GUNUNG_SEGARA,
        "gunung_segara",
        towgs84!(-403., 684., 41.),
        BESSEL,
        "Gunung Segara Jakarta",
    );
    datum!(
        RNB72,
        "rnb72",
        towgs84!(106.869, -52.2978, 103.724, -0.33657, 0.456955, -1.84218, 1.),
        INTL,
        "Reseau National Belge 1972",
    );

    /// Static datums table
    pub(super) const DATUMS: [&DatumDefn; 17] = [
        &WGS84,
        &GGRS87,
        &NAD83,
        &NAD27,
        &POTSDAM,
        &CARTHAGE,
        &HERMANNSKOGEL,
        &IRE65,
        &NZGD49,
        &OSGB36,
        &CH1903,
        &OSNI52,
        &RASSADIRAN,
        &S_JTSK,
        &BEDUARAM,
        &GUNUNG_SEGARA,
        &RNB72,
    ];
}

/// Return the datum definition
pub fn find_datum(name: &str) -> Option<&DatumDefn> {
    constants::DATUMS
        .iter()
        .find(|d| d.id.eq_ignore_ascii_case(name))
        .copied()
}
