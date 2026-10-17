//!
//! Wasm bindgen entry point
//!
//! All items are exported in the `geodesic` js namespace.
//!
use crate::geodesic::{Flags, Mask};
use wasm_bindgen::prelude::*;

// ----------------------------
// Masks and flags
// ---------------------------

/// Bit masks for the quantities to compute (or the capabilities of a line).
///
/// Values may be combined with `|`.
#[wasm_bindgen(js_namespace = geodesic, js_name = Mask)]
#[derive(Clone, Copy)]
pub enum JsMask {
    None = 0,
    Latitude = 128,
    Longitude = 264,
    Azimuth = 512,
    Distance = 1025,
    DistanceIn = 2051,
    ReducedLength = 4101,
    GeodesicScale = 8197,
    Area = 16400,
    All = 32671,
}

/// Flags for the direct problem and positions along a line.
///
/// Values may be combined with `|`.
#[wasm_bindgen(js_namespace = geodesic, js_name = Flags)]
#[derive(Clone, Copy)]
pub enum JsFlags {
    None = 0,
    ArcMode = 1,
    LongUnroll = 32768,
}

// Keep js values in sync with Rust masks
const _: () = {
    assert!(JsMask::None as u32 == Mask::NONE.0);
    assert!(JsMask::Latitude as u32 == Mask::LATITUDE.0);
    assert!(JsMask::Longitude as u32 == Mask::LONGITUDE.0);
    assert!(JsMask::Azimuth as u32 == Mask::AZIMUTH.0);
    assert!(JsMask::Distance as u32 == Mask::DISTANCE.0);
    assert!(JsMask::DistanceIn as u32 == Mask::DISTANCE_IN.0);
    assert!(JsMask::ReducedLength as u32 == Mask::REDUCED_LENGTH.0);
    assert!(JsMask::GeodesicScale as u32 == Mask::GEODESIC_SCALE.0);
    assert!(JsMask::Area as u32 == Mask::AREA.0);
    assert!(JsMask::All as u32 == Mask::ALL.0);
};

fn to_mask(bits: u32) -> Mask {
    Mask(bits & Mask::ALL.0)
}

fn to_flags(bits: u32) -> Flags {
    let mut flags = Flags::NONE;
    if bits & JsFlags::ArcMode as u32 != 0 {
        flags = flags | Flags::ARC_MODE;
    }
    if bits & JsFlags::LongUnroll as u32 != 0 {
        flags = flags | Flags::LONG_UNROLL;
    }
    flags
}

// ----------------------------
// Results
// ---------------------------

/// Result of the direct problem: `[lat2, lon2, azi2]`
#[wasm_bindgen(js_namespace = geodesic)]
#[derive(Clone, Copy)]
pub struct Direct {
    pub lat2: f64,
    pub lon2: f64,
    pub azi2: f64,
}

/// Result of the inverse problem
#[wasm_bindgen(js_namespace = geodesic)]
#[derive(Clone, Copy)]
pub struct Inverse {
    /// Distance between the points [meters]
    pub s12: f64,
    /// Azimuth at 1st point [degrees]
    pub azi1: f64,
    /// Azimuth at 2nd point [degrees]
    pub azi2: f64,
}

/// Result of the general inverse problem
///
/// Quantities not requested in the mask are set to `NaN`.
#[wasm_bindgen(js_namespace = geodesic)]
#[derive(Clone, Copy)]
pub struct GenInverse {
    /// Arc length between the points [degrees]
    pub a12: f64,
    /// Distance between the points [meters]
    pub s12: f64,
    /// Azimuth at 1st point [degrees]
    pub azi1: f64,
    /// Azimuth at 2nd point [degrees]
    pub azi2: f64,
    /// Reduced length [meters]
    pub m12: f64,
    /// Geodesic scale of point 2 relative to point 1
    pub scale12: f64,
    /// Geodesic scale of point 1 relative to point 2
    pub scale21: f64,
    /// Area under the geodesic [meters²]
    pub area12: f64,
}

impl From<crate::Inverse> for GenInverse {
    fn from(r: crate::Inverse) -> Self {
        Self {
            a12: r.a12,
            s12: r.s12,
            azi1: r.azi1,
            azi2: r.azi2,
            m12: r.m12,
            scale12: r.scale12,
            scale21: r.scale21,
            area12: r.area12,
        }
    }
}

/// Position on a geodesic
///
/// Quantities not requested in the mask are set to `NaN`.
#[wasm_bindgen(js_namespace = geodesic)]
#[derive(Clone, Copy)]
pub struct Position {
    /// Arc length from 1st point [degrees]
    pub a12: f64,
    /// Latitude of 2nd point [degrees]
    pub lat2: f64,
    /// Longitude of 2nd point [degrees]
    pub lon2: f64,
    /// Azimuth at 2nd point [degrees]
    pub azi2: f64,
    /// Distance from 1st point [meters]
    pub s12: f64,
    /// Reduced length [meters]
    pub m12: f64,
    /// Geodesic scale of point 2 relative to point 1
    pub scale12: f64,
    /// Geodesic scale of point 1 relative to point 2
    pub scale21: f64,
    /// Area under the geodesic [meters²]
    pub area12: f64,
}

impl From<crate::Position> for Position {
    fn from(p: crate::Position) -> Self {
        Self {
            a12: p.a12,
            lat2: p.lat2,
            lon2: p.lon2,
            azi2: p.azi2,
            s12: p.s12,
            m12: p.m12,
            scale12: p.scale12,
            scale21: p.scale21,
            area12: p.area12,
        }
    }
}

/// Result of polygon computations
#[wasm_bindgen(js_namespace = geodesic)]
#[derive(Clone, Copy)]
pub struct PolygonResult {
    /// The number of points
    pub num: u32,
    /// The perimeter of the polygon or the length of the polyline [meters]
    pub perimeter: f64,
    /// The area of the polygon [meters²], `NaN` for a polyline
    pub area: f64,
}

impl From<crate::PolygonResult> for PolygonResult {
    fn from(r: crate::PolygonResult) -> Self {
        Self {
            num: r.num,
            perimeter: r.perimeter,
            area: r.area.unwrap_or(f64::NAN),
        }
    }
}

// ----------------------------
// Wrapper for Geodesic
// ---------------------------

/// Ellipsoid on which geodesic calculations are computed
#[wasm_bindgen(js_namespace = geodesic)]
pub struct Geodesic {
    inner: crate::Geodesic,
}

#[wasm_bindgen(js_namespace = geodesic)]
impl Geodesic {
    /// Create new ellipsoid with semi-major axis `a` in meters and a flattening `f`
    #[wasm_bindgen(constructor)]
    pub fn new(a: f64, f: f64) -> Geodesic {
        Self {
            inner: crate::Geodesic::new(a, f),
        }
    }

    /// Create the WGS84 ellipsoid
    pub fn wgs84() -> Geodesic {
        Self {
            inner: crate::Geodesic::wgs84(),
        }
    }

    /// The equatorial radius [meters]
    #[wasm_bindgen(getter)]
    pub fn a(&self) -> f64 {
        self.inner.a()
    }

    /// The flattening
    #[wasm_bindgen(getter)]
    pub fn f(&self) -> f64 {
        self.inner.f()
    }

    /// Solve the direct geodesic problem
    pub fn direct(&self, lat1: f64, lon1: f64, azi1: f64, s12: f64) -> Direct {
        let (lat2, lon2, azi2) = self.inner.direct(lat1, lon1, azi1, s12);
        Direct { lat2, lon2, azi2 }
    }

    /// The general direct geodesic problem
    #[wasm_bindgen(js_name = genDirect)]
    pub fn gen_direct(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        flags: u32,
        s12_a12: f64,
        outmask: u32,
    ) -> Position {
        self.inner
            .gen_direct(lat1, lon1, azi1, to_flags(flags), s12_a12, to_mask(outmask))
            .into()
    }

    /// Solve the inverse geodesic problem
    pub fn inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> Inverse {
        let (s12, azi1, azi2) = self.inner.inverse(lat1, lon1, lat2, lon2);
        Inverse { s12, azi1, azi2 }
    }

    /// The general inverse geodesic problem
    #[wasm_bindgen(js_name = genInverse)]
    pub fn gen_inverse(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        outmask: u32,
    ) -> GenInverse {
        self.inner
            .gen_inverse(lat1, lon1, lat2, lon2, to_mask(outmask))
            .into()
    }

    /// Create a geodesic line starting at (`lat1`, `lon1`) with azimuth `azi1`
    pub fn line(&self, lat1: f64, lon1: f64, azi1: f64, caps: u32) -> GeodesicLine {
        GeodesicLine {
            inner: self.inner.line(lat1, lon1, azi1, to_mask(caps)),
        }
    }

    /// Create a geodesic line in terms of the direct geodesic problem
    #[wasm_bindgen(js_name = directLine)]
    pub fn direct_line(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        s12: f64,
        caps: u32,
    ) -> GeodesicLine {
        GeodesicLine {
            inner: self.inner.direct_line(lat1, lon1, azi1, s12, to_mask(caps)),
        }
    }

    /// Create a geodesic line in terms of the general direct geodesic problem
    #[wasm_bindgen(js_name = genDirectLine)]
    pub fn gen_direct_line(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        flags: u32,
        s12_a12: f64,
        caps: u32,
    ) -> GeodesicLine {
        GeodesicLine {
            inner: self.inner.gen_direct_line(
                lat1,
                lon1,
                azi1,
                to_flags(flags),
                s12_a12,
                to_mask(caps),
            ),
        }
    }

    /// Create a geodesic line in terms of the inverse geodesic problem
    #[wasm_bindgen(js_name = inverseLine)]
    pub fn inverse_line(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        caps: u32,
    ) -> GeodesicLine {
        GeodesicLine {
            inner: self
                .inner
                .inverse_line(lat1, lon1, lat2, lon2, to_mask(caps)),
        }
    }

    /// Create a new polygon (or polyline if `polyline` is true)
    pub fn polygon(&self, polyline: bool) -> Polygon {
        Polygon::new(self, polyline)
    }

    /// Compute the area and perimeter of a geodesic polygon
    ///
    /// `lats` and `lons` are the coordinates of the vertices [degrees].
    /// The area is signed with counter-clockwise traversal being
    /// treated as positive.
    #[wasm_bindgen(js_name = polygonArea)]
    pub fn polygon_area(&self, lats: &[f64], lons: &[f64]) -> Result<PolygonResult, JsError> {
        if lats.len() != lons.len() {
            return Err(JsError::new("lats and lons must have the same length"));
        }
        let mut p = self.inner.polygon(false);
        lats.iter()
            .zip(lons)
            .for_each(|(lat, lon)| p.add_point(*lat, *lon));
        Ok(p.compute(false, true).into())
    }
}

// ----------------------------
// Wrapper for GeodesicLine
// ---------------------------

/// A geodesic line, for efficient computation of many points
/// along a single geodesic
#[wasm_bindgen(js_namespace = geodesic)]
pub struct GeodesicLine {
    inner: crate::GeodesicLine,
}

#[wasm_bindgen(js_namespace = geodesic)]
impl GeodesicLine {
    /// Latitude of point 1 [degrees]
    #[wasm_bindgen(getter)]
    pub fn lat1(&self) -> f64 {
        self.inner.lat1()
    }

    /// Longitude of point 1 [degrees]
    #[wasm_bindgen(getter)]
    pub fn lon1(&self) -> f64 {
        self.inner.lon1()
    }

    /// Azimuth at point 1 [degrees]
    #[wasm_bindgen(getter)]
    pub fn azi1(&self) -> f64 {
        self.inner.azi1()
    }

    /// Arc length to point 3 [degrees]
    #[wasm_bindgen(getter)]
    pub fn a13(&self) -> f64 {
        self.inner.a13()
    }

    /// Distance to point 3 [meters]
    #[wasm_bindgen(getter)]
    pub fn s13(&self) -> f64 {
        self.inner.s13()
    }

    /// The capabilities of the line
    #[wasm_bindgen(getter)]
    pub fn caps(&self) -> u32 {
        self.inner.caps().0
    }

    /// Compute the position at distance `s12` [meters] from point 1
    pub fn position(&self, s12: f64) -> Direct {
        let (lat2, lon2, azi2) = self.inner.position(s12);
        Direct { lat2, lon2, azi2 }
    }

    /// The general position function
    #[wasm_bindgen(js_name = genPosition)]
    pub fn gen_position(&self, flags: u32, s12_a12: f64, outmask: u32) -> Position {
        self.inner
            .gen_position(to_flags(flags), s12_a12, to_mask(outmask))
            .into()
    }

    /// Set the distance to point 3 [meters]
    #[wasm_bindgen(js_name = setDistance)]
    pub fn set_distance(&mut self, s13: f64) {
        self.inner.set_distance(s13)
    }

    /// Set the arc length to point 3 [degrees]
    #[wasm_bindgen(js_name = setArc)]
    pub fn set_arc(&mut self, a13: f64) {
        self.inner.set_arc(a13)
    }
}

// ----------------------------
// Wrapper for Polygon
// ---------------------------

/// Accumulate the perimeter and area of a geodesic polygon (or the
/// length of a polyline).
#[wasm_bindgen(js_namespace = geodesic)]
pub struct Polygon {
    inner: crate::Polygon<crate::Geodesic>,
}

#[wasm_bindgen(js_namespace = geodesic)]
impl Polygon {
    /// Create a new polygon on the ellipsoid `g`
    #[wasm_bindgen(constructor)]
    pub fn new(g: &Geodesic, polyline: bool) -> Polygon {
        Self {
            inner: crate::Polygon::new(g.inner.clone(), polyline),
        }
    }

    /// Clear the polygon, allowing a new polygon to be started
    pub fn clear(&mut self) {
        self.inner.clear()
    }

    /// The number of points so far
    #[wasm_bindgen(getter)]
    pub fn num(&self) -> u32 {
        self.inner.num()
    }

    /// Add a point (`lat`, `lon`) [degrees]
    #[wasm_bindgen(js_name = addPoint)]
    pub fn add_point(&mut self, lat: f64, lon: f64) {
        self.inner.add_point(lat, lon)
    }

    /// Add an edge given by the azimuth `azi` [degrees] at the current
    /// point and the distance `s` [meters] to the next point
    ///
    /// This does nothing if no points have been added yet.
    #[wasm_bindgen(js_name = addEdge)]
    pub fn add_edge(&mut self, azi: f64, s: f64) {
        self.inner.add_edge(azi, s)
    }

    /// Return the results for the polygon
    pub fn compute(&self, reverse: bool, sign: bool) -> PolygonResult {
        self.inner.compute(reverse, sign).into()
    }

    /// Return the results assuming a tentative final test point is added
    #[wasm_bindgen(js_name = testPoint)]
    pub fn test_point(&self, lat: f64, lon: f64, reverse: bool, sign: bool) -> PolygonResult {
        self.inner.test_point(lat, lon, reverse, sign).into()
    }

    /// Return the results assuming a tentative final test edge is added
    #[wasm_bindgen(js_name = testEdge)]
    pub fn test_edge(&self, azi: f64, s: f64, reverse: bool, sign: bool) -> PolygonResult {
        self.inner.test_edge(azi, s, reverse, sign).into()
    }
}
