//!
//! Proj4 Ellipsoids definitions
//!

/// Ellipsoid flattening may be defined either by
/// the knowledge of its minor axis or by its reverse flattening
pub enum FlatteningParam {
    MinorAxis(f64),
    InvFlat(f64),
}

use FlatteningParam::*;

pub struct EllipsoidDefn {
    pub id: &'static str,
    pub a: f64,
    pub rf_or_b: FlatteningParam,
    //pub comment: &'static str,
}

#[rustfmt::skip]
pub mod constants {
    use super::*;

    macro_rules! ellps {
        ($name:ident, $id:expr, $a:expr, $f:expr, $c:expr) => {
            pub(crate) const $name: EllipsoidDefn = EllipsoidDefn {
                id: $id,
                a: $a,
                rf_or_b: $f,
                //comment: $c,
            };
        };
    }

    ellps!(MERIT,     "MERIT",     6_378_137.,     InvFlat(298.257),           "MERIT 1983");
    ellps!(SGS85,     "SGS85",     6_378_136.,     InvFlat(298.257),           "Soviet Geodetic System 85");
    ellps!(GRS80,     "GRS80",     6_378_137.,     InvFlat(298.257222101),     "GRS 1980(IUGG, 1980)");
    ellps!(IAU76,     "IAU76",     6_378_140.,     InvFlat(298.257),           "IAU 1976");
    ellps!(AIRY,      "airy",      6_377_563.396,  InvFlat(299.3249646),       "Airy 1830");
    ellps!(APL4_9,    "APL4.9",    6_378_137.,     InvFlat(298.25),            "Appl. Physics. 1965");
    ellps!(NWL9D,     "NWL9D",     6_378_145.,     InvFlat(298.25),            "Naval Weapons Lab., 1965");
    ellps!(MOD_AIRY,  "mod_airy",  6_377_340.189,  MinorAxis(6_356_034.446),   "Modified Airy");
    ellps!(ANDRAE,    "andrae",    6_377_104.43,   InvFlat(300.0),             "Andrae 1876 (Den., Iclnd.)");
    ellps!(DANISH,    "danish",    6_377_019.256_3, InvFlat(300.0),             "Andrae 1876 (Denmark, Iceland)");
    ellps!(AUST_SA,   "aust_SA",   6_378_160.,     InvFlat(298.25),            "Australian Natl & S. Amer. 1969");
    ellps!(GRS67,     "GRS67",     6_378_160.,     InvFlat(298.2471674270),    "GRS 67(IUGG 1967)");
    ellps!(GSK2011,   "GSK2011",   6_378_136.5,    InvFlat(298.2564151),       "GSK-2011");
    ellps!(BESSEL,    "bessel",    6_377_397.155,  InvFlat(299.1528128),       "Bessel 1841");
    ellps!(BESS_NAM,  "bess_nam",  6_377_483.865,  InvFlat(299.1528128),       "Bessel 1841 (Namibia)");
    ellps!(CLRK66,    "clrk66",    6_378_206.4,    MinorAxis(6_356_583.8),     "Clarke 1866");
    ellps!(CLRK80,    "clrk80",    6_378_249.145,  InvFlat(293.4663),          "Clarke 1880 mod.");
    ellps!(CLRK80IGN, "clrk80ign", 6_378_249.2,    InvFlat(293.4660212936269), "Clarke 1880 (IGN).");
    ellps!(CPM,       "CPM",       6_375_738.7,    InvFlat(334.29),            "Comm. des Poids et Mesures 1799");
    ellps!(DELMBR,    "delmbr",    6_376_428.,     InvFlat(311.5),             "Delambre 1810 (Belgium)");
    ellps!(ENGELIS,   "engelis",   6_378_136.05,   InvFlat(298.2566),          "Engelis 1985");
    ellps!(EVRST30,   "evrst30",   6_377_276.345,  InvFlat(300.8017),          "Everest 1830");
    ellps!(EVRST48,   "evrst48",   6_377_304.063,  InvFlat(300.8017),          "Everest 1948");
    ellps!(EVRST56,   "evrst56",   6_377_301.243,  InvFlat(300.8017),          "Everest 1956");
    ellps!(EVRST69,   "evrst69",   6_377_295.664,  InvFlat(300.8017),          "Everest 1969");
    ellps!(EVRSTSS,   "evrstSS",   6_377_298.556,  InvFlat(300.8017),          "Everest (Sabah & Sarawak)");
    ellps!(FSCHR60,   "fschr60",   6_378_166.,     InvFlat(298.3),             "Fischer (Mercury Datum) 1960");
    ellps!(FSCHR60M,  "fschr60m",  6_378_155.,     InvFlat(298.3),             "Modified Fischer 1960");
    ellps!(FSCHR68,   "fschr68",   6_378_150.,     InvFlat(298.3),             "Fischer 1968");
    ellps!(HELMERT,   "helmert",   6_378_200.,     InvFlat(298.3),             "Helmert 1906");
    ellps!(HOUGH,     "hough",     6_378_270.,     InvFlat(297.),              "Hough");
    ellps!(INTL,      "intl",      6_378_388.,     InvFlat(297.),              "International 1924 (Hayford 1909, 1910)");
    ellps!(KRASS,     "krass",     6_378_245.,     InvFlat(298.3),             "Krassovsky, 1942");
    ellps!(KAULA,     "kaula",     6_378_163.,     InvFlat(298.24),            "Kaula 1961");
    ellps!(LERCH,     "lerch",     6_378_139.,     InvFlat(298.257),           "Lerch 1979");
    ellps!(MPRTS,     "mprts",     6_397_300.,     InvFlat(191.),              "Maupertius 1738");
    ellps!(NEW_INTL,  "new_intl",  6_378_157.5,    MinorAxis(6_356_772.2),     "New International 1967");
    ellps!(PLESSIS,   "plessis",   6_376_523.,     MinorAxis(6_355_863.),      "Plessis 1817 (France)");
    ellps!(PZ90,      "PZ90",      6_378_136.,     InvFlat(298.25784),          "PZ-90");
    ellps!(SEASIA,    "SEasia",    6_378_155.,     MinorAxis(6_356_773.320_5),  "Southeast Asia");
    ellps!(WALBECK,   "walbeck",   6_376_896.,     MinorAxis(6_355_834.846_7),  "Walbeck");
    ellps!(WGS60,     "WGS60",     6_378_165.,     InvFlat(298.3),             "WGS 60");
    ellps!(WGS66,     "WGS66",     6_378_145.,     InvFlat(298.25),            "WGS 66");
    ellps!(WGS72,     "WGS72",     6_378_135.,     InvFlat(298.26),            "WGS 72");
    ellps!(WGS84,     "WGS84",     6_378_137.,     InvFlat(298.257_223_563),   "WGS 84");
    ellps!(SPHERE,    "sphere",    6_370_997.,     MinorAxis(6_370_997.),      "Normal Sphere (r=6370997)");

    ///
    /// Static ellipsoids table
    ///
    /// Format: (id, major axis (a), FlatteningParam (b or rf), comment)
    pub (super) const ELLIPSOIDS: [&EllipsoidDefn;46] = [
        &MERIT,
        &SGS85,
        &GRS80,
        &IAU76,
        &AIRY,
        &APL4_9,
        &NWL9D,
        &MOD_AIRY,
        &ANDRAE,
        &DANISH,
        &AUST_SA,
        &GRS67,
        &GSK2011,
        &BESSEL,
        &BESS_NAM,
        &CLRK66,
        &CLRK80,
        &CLRK80IGN,
        &CPM,
        &DELMBR,
        &ENGELIS,
        &EVRST30,
        &EVRST48,
        &EVRST56,
        &EVRST69,
        &EVRSTSS,
        &FSCHR60,
        &FSCHR60M,
        &FSCHR68,
        &HELMERT,
        &HOUGH,
        &INTL,
        &KRASS,
        &KAULA,
        &LERCH,
        &MPRTS,
        &NEW_INTL,
        &PLESSIS,
        &PZ90,
        &SEASIA,
        &WALBECK,
        &WGS60,
        &WGS66,
        &WGS72,
        &WGS84,
        &SPHERE,
    ];

}

/// Return the ellipse definition
pub fn find_ellipsoid(name: &str) -> Option<&EllipsoidDefn> {
    constants::ELLIPSOIDS
        .iter()
        .find(|e| e.id.eq_ignore_ascii_case(name))
        .copied()
}
