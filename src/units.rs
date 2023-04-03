//!
//! Predefined units for conversion
//!

#[derive(Debug, Copy, Clone)]
pub struct UnitDefn {
    pub name: &'static str,
    pub to_meter: f64,
}

macro_rules! unit {
    ($name:expr, $display:expr, $comment:expr, $to_meter:expr) => {
        UnitDefn {
            name: $name,
            to_meter: $to_meter,
        }
    };
}

pub const METER: UnitDefn = unit!("m", "1", "Meter", 1.0);

pub const DEGREES: &str = "degrees";

mod constants {
    use super::*;
    /// Static units table
    /// id, to_meter, display to_meter value, comment, to_meter
    #[rustfmt::skip]
    pub const UNITS: [UnitDefn;21] = [
        unit!("km",      "1000",                 "Kilometer",                    1000.0),
        unit!("m",       "1",                    "Meter",                        1.0),
        unit!("dm",      "1/10",                 "Decimeter",                    0.1),
        unit!("cm",      "1/100",                "Centimeter",                   0.01),
        unit!("mm",      "1/1000",               "Millimeter",                   0.001),
        unit!("kmi",     "1852",                 "International Nautical Mile",  1852.0),
        unit!("in",      "0.0254",               "International Inch",           0.0254),
        unit!("ft",      "0.3048",               "International Foot",           0.3048),
        unit!("yd",      "0.9144",               "International Yard",           0.9144),
        unit!("mi",      "1609.344",             "International Statute Mile",   1609.344),
        unit!("fath",    "1.8288",               "International Fathom",         1.8288),
        unit!("ch",      "20.1168",              "International Chain",          20.1168),
        unit!("link",    "0.201168",             "International Link",           0.201168),
        unit!("us-in",   "1/39.37",              "U.S. Surveyor's Inch",         100./3937.0),
        unit!("us-ft",   "0.304800609601219",    "U.S. Surveyor's Foot",         1200./3937.0),
        unit!("us-yd",   "0.914401828803658",    "U.S. Surveyor's Yard",         3600./3937.0),
        unit!("us-ch",   "20.11684023368047",    "U.S. Surveyor's Chain",        79200./3937.0),
        unit!("us-mi",   "1609.347218694437",    "U.S. Surveyor's Statute Mile", 6336000./3937.0),
        unit!("ind-yd",  "0.91439523",           "Indian Yard",                  0.91439523),
        unit!("ind-ft",  "0.30479841",           "Indian Foot",                  0.30479841),
        unit!("ind-ch",  "20.11669506",          "Indian Chain",                 20.11669506),
    ];
}

pub fn from_value(to_meter: f64) -> UnitDefn {
    UnitDefn { name: "", to_meter }
}

/// Return the unit definition
pub fn find_units(name: &str) -> Option<UnitDefn> {
    constants::UNITS
        .iter()
        .find(|d| d.name.eq_ignore_ascii_case(name))
        .copied()
}
