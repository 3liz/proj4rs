//! Perimeter and area of geodesic polygons
use crate::geodesic::{Flags, Geodesic, Mask};
use crate::math::*;
use std::borrow::Borrow;
use std::f64::consts::PI;

/// Result of polygon computations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolygonResult {
    /// The number of points
    pub num: u32,
    /// The perimeter of the polygon or the length of the polyline \[meters\]
    pub perimeter: f64,
    /// The area of the polygon \[meters²\], `None` for a polyline
    pub area: Option<f64>,
}

/// Accumulate the perimeter and area of a geodesic polygon (or the length
/// of a polyline).
///
/// The area and perimeter are accumulated at two times the standard floating
/// point precision to guard against the loss of accuracy with many-sided
/// polygons.  At any point you can ask for the perimeter and area so far.
///
/// Arbitrarily complex polygons are allowed.  In the case of
/// self-intersecting polygons the area is accumulated "algebraically", e.g.,
/// the areas of the 2 loops in a figure-8 polygon will partially cancel.
/// There's no need to "close" the polygon by repeating the first vertex.
///
/// ```rust
/// // Compute the perimeter and area of the geodesic triangle with
/// // vertices (0°N,0°E), (0°N,90°E), (90°N,0°E).
/// use proj4rs_geodesic::Geodesic;
/// let g = Geodesic::wgs84();
/// let mut p = g.polygon(false);
/// p.add_point(0., 0.);
/// p.add_point(0., 90.);
/// p.add_point(90., 0.);
/// let r = p.compute(false, true);
/// assert_eq!(r.num, 3);
/// assert_eq!(r.perimeter, 30022685.630020067);
/// assert_eq!(r.area, Some(63758202715511.055));
/// ```
///
/// The ellipsoid is held through any type implementing
/// [`Borrow<Geodesic>`](Borrow): a reference (as returned by
/// [`Geodesic::polygon`]), an owned [`Geodesic`] or a shared pointer
/// such as [`Rc<Geodesic>`](std::rc::Rc).
///
/// ```rust
/// use proj4rs_geodesic::{Geodesic, Polygon};
/// use std::rc::Rc;
/// let g = Rc::new(Geodesic::wgs84());
/// let mut p = Polygon::new(Rc::clone(&g), false);
/// p.add_point(0., 0.);
/// p.add_point(0., 90.);
/// p.add_point(90., 0.);
/// assert_eq!(p.compute(false, true).area, Some(63758202715511.055));
/// ```
#[derive(Debug, Clone)]
pub struct Polygon<G: Borrow<Geodesic>> {
    g: G,
    lat: f64,
    lon: f64,
    lat0: f64,
    lon0: f64,
    area: Accumulator,
    perimeter: Accumulator,
    polyline: bool,
    crossings: i32,
    num: u32,
}

impl<G: Borrow<Geodesic>> Polygon<G> {
    /// Create a new polygon on the ellipsoid `g`.
    ///
    /// If `polyline` is true, then the vertices and edges define a polyline
    /// and only the perimeter is computed.
    pub fn new(g: G, polyline: bool) -> Self {
        Self {
            g,
            lat: f64::NAN,
            lon: f64::NAN,
            lat0: f64::NAN,
            lon0: f64::NAN,
            area: Accumulator::default(),
            perimeter: Accumulator::default(),
            polyline,
            crossings: 0,
            num: 0,
        }
    }

    /// Clear the polygon, allowing a new polygon to be started.
    pub fn clear(&mut self) {
        self.lat = f64::NAN;
        self.lon = f64::NAN;
        self.lat0 = f64::NAN;
        self.lon0 = f64::NAN;
        self.area = Accumulator::default();
        self.perimeter = Accumulator::default();
        self.crossings = 0;
        self.num = 0;
    }

    /// The ellipsoid of the polygon
    #[inline(always)]
    pub fn geodesic(&self) -> &Geodesic {
        self.g.borrow()
    }

    /// Return true if this is a polyline
    #[inline(always)]
    pub fn is_polyline(&self) -> bool {
        self.polyline
    }

    /// The number of points so far
    #[inline(always)]
    pub fn num(&self) -> u32 {
        self.num
    }

    /// The current point (`lat`, `lon`) \[degrees\]
    #[inline(always)]
    pub fn current_point(&self) -> (f64, f64) {
        (self.lat, self.lon)
    }

    fn outmask(&self) -> Mask {
        if self.polyline {
            Mask::DISTANCE
        } else {
            Mask::DISTANCE | Mask::AREA
        }
    }

    /// Add a point (`lat`, `lon`) \[degrees\] to the polygon or polyline.
    pub fn add_point(&mut self, lat: f64, lon: f64) {
        if self.num == 0 {
            self.lat0 = lat;
            self.lon0 = lon;
        } else {
            let r = self
                .geodesic()
                .gen_inverse(self.lat, self.lon, lat, lon, self.outmask());
            self.perimeter.add(r.s12);
            if !self.polyline {
                self.area.add(r.area12);
                self.crossings += transit(self.lon, lon);
            }
        }
        self.lat = lat;
        self.lon = lon;
        self.num += 1;
    }

    /// Add an edge to the polygon or polyline, given by the azimuth
    /// `azi` \[degrees\] at the current point and the distance `s` \[meters\]
    /// to the next point.
    ///
    /// This does nothing if no points have been added yet.
    pub fn add_edge(&mut self, azi: f64, s: f64) {
        if self.num > 0 {
            let outmask = if self.polyline {
                Mask::LATITUDE | Mask::LONGITUDE
            } else {
                Mask::LATITUDE | Mask::LONGITUDE | Mask::AREA
            };
            let p =
                self.geodesic()
                    .gen_direct(self.lat, self.lon, azi, Flags::LONG_UNROLL, s, outmask);
            self.perimeter.add(s);
            if !self.polyline {
                self.area.add(p.area12);
                self.crossings += transitdirect(self.lon, p.lon2);
            }
            self.lat = p.lat2;
            self.lon = p.lon2;
            self.num += 1;
        }
    }

    fn area0(&self) -> f64 {
        4.0 * PI * self.geodesic().c2
    }

    /// Return the results for the polygon.
    ///
    /// If `reverse` is true then clockwise (instead of counter-clockwise)
    /// traversal counts as a positive area. If `sign` is true then return a
    /// signed result for the area if the polygon is traversed in the "wrong"
    /// direction instead of returning the area for the rest of the earth.
    ///
    /// More points can be added to the polygon after this call.
    pub fn compute(&self, reverse: bool, sign: bool) -> PolygonResult {
        let num = self.num;
        if num < 2 {
            return PolygonResult {
                num,
                perimeter: 0.0,
                area: (!self.polyline).then_some(0.0),
            };
        }
        if self.polyline {
            return PolygonResult {
                num,
                perimeter: self.perimeter.sum(),
                area: None,
            };
        }
        let r = self.geodesic().gen_inverse(
            self.lat,
            self.lon,
            self.lat0,
            self.lon0,
            Mask::DISTANCE | Mask::AREA,
        );
        let mut area = self.area;
        area.add(r.area12);
        PolygonResult {
            num,
            perimeter: self.perimeter.sum_with(r.s12),
            area: Some(area_reduce_acc(
                area,
                self.area0(),
                self.crossings + transit(self.lon, self.lon0),
                reverse,
                sign,
            )),
        }
    }

    /// Return the results assuming a tentative final test point is added;
    /// however, the data for the test point is not saved.
    ///
    /// This lets you report a running result for the perimeter and area as
    /// the user moves the mouse cursor.  Ordinary floating point arithmetic is
    /// used to accumulate the data for the test point; thus the area and
    /// perimeter returned are less accurate than if
    /// [`Polygon::add_point`] and [`Polygon::compute`] are used.
    pub fn test_point(&self, lat: f64, lon: f64, reverse: bool, sign: bool) -> PolygonResult {
        let num = self.num + 1;
        if num == 1 {
            return PolygonResult {
                num,
                perimeter: 0.0,
                area: (!self.polyline).then_some(0.0),
            };
        }
        let mut perimeter = self.perimeter.sum();
        let mut tempsum = if self.polyline { 0.0 } else { self.area.sum() };
        let mut crossings = self.crossings;
        let edges: &[_] = if self.polyline {
            &[(self.lat, self.lon, lat, lon)]
        } else {
            &[
                (self.lat, self.lon, lat, lon),
                (lat, lon, self.lat0, self.lon0),
            ]
        };
        for &(lat1, lon1, lat2, lon2) in edges {
            let r = self
                .geodesic()
                .gen_inverse(lat1, lon1, lat2, lon2, self.outmask());
            perimeter += r.s12;
            if !self.polyline {
                tempsum += r.area12;
                crossings += transit(lon1, lon2);
            }
        }

        PolygonResult {
            num,
            perimeter,
            area: (!self.polyline)
                .then(|| area_reduce(tempsum, self.area0(), crossings, reverse, sign)),
        }
    }

    /// Return the results assuming a tentative final test point is added via
    /// an azimuth `azi` \[degrees\] and distance `s` \[meters\]; however, the
    /// data for the test point is not saved.
    ///
    /// This lets you report a running result for the perimeter and area as
    /// the user moves the mouse cursor.  Ordinary floating point arithmetic is
    /// used to accumulate the data for the test point; thus the area and
    /// perimeter returned are less accurate than if
    /// [`Polygon::add_edge`] and [`Polygon::compute`] are used.
    pub fn test_edge(&self, azi: f64, s: f64, reverse: bool, sign: bool) -> PolygonResult {
        let num = self.num + 1;
        if num == 1 {
            // we don't have a starting point!
            return PolygonResult {
                num: 0,
                perimeter: f64::NAN,
                area: (!self.polyline).then_some(f64::NAN),
            };
        }
        let mut perimeter = self.perimeter.sum() + s;
        if self.polyline {
            return PolygonResult {
                num,
                perimeter,
                area: None,
            };
        }

        let mut tempsum = self.area.sum();
        let mut crossings = self.crossings;

        let p = self.geodesic().gen_direct(
            self.lat,
            self.lon,
            azi,
            Flags::LONG_UNROLL,
            s,
            Mask::LATITUDE | Mask::LONGITUDE | Mask::AREA,
        );
        tempsum += p.area12;
        crossings += transitdirect(self.lon, p.lon2);
        let r = self.geodesic().gen_inverse(
            p.lat2,
            p.lon2,
            self.lat0,
            self.lon0,
            Mask::DISTANCE | Mask::AREA,
        );
        perimeter += r.s12;
        tempsum += r.area12;
        crossings += transit(p.lon2, self.lon0);

        PolygonResult {
            num,
            perimeter,
            area: Some(area_reduce(tempsum, self.area0(), crossings, reverse, sign)),
        }
    }
}

impl Geodesic {
    /// Create a new [`Polygon`] (or polyline if `polyline` is true).
    pub fn polygon(&self, polyline: bool) -> Polygon<&Geodesic> {
        Polygon::new(self, polyline)
    }

    /// A simple interface for computing the area of a geodesic polygon.
    ///
    /// `points` are the (`lat`, `lon`) vertices of the polygon \[degrees\].
    /// Returns the `(area, perimeter)` of the polygon in meters² and meters.
    ///
    /// Arbitrarily complex polygons are allowed.  In the case self-intersecting
    /// of polygons the area is accumulated "algebraically", e.g., the areas of
    /// the 2 loops in a figure-8 polygon will partially cancel.  There's no need
    /// to "close" the polygon by repeating the first vertex.  The area returned
    /// is signed with counter-clockwise traversal being treated as positive.
    ///
    /// ```rust
    /// // Compute the area of Antarctica
    /// use proj4rs_geodesic::Geodesic;
    /// let lats = [-72.9, -71.9, -74.9, -74.3, -77.5, -77.4, -71.7, -65.9, -65.7,
    ///             -66.6, -66.9, -69.8, -70.0, -71.0, -77.3, -77.9, -74.7];
    /// let lons = [-74., -102., -102., -131., -163., 163., 172., 140., 113.,
    ///             88., 59., 25., -4., -14., -33., -46., -61.];
    /// let g = Geodesic::wgs84();
    /// let (area, perimeter) = g.polygon_area(lats.into_iter().zip(lons));
    /// assert_eq!(area, 13376856682207.375);
    /// assert_eq!(perimeter, 14710425.40697364);
    /// ```
    pub fn polygon_area<I>(&self, points: I) -> (f64, f64)
    where
        I: IntoIterator<Item = (f64, f64)>,
    {
        let mut p = self.polygon(false);
        for (lat, lon) in points {
            p.add_point(lat, lon);
        }
        let r = p.compute(false, true);
        (r.area.unwrap_or(0.0), r.perimeter)
    }
}

/// Return 1 or -1 if crossing prime meridian in east or west direction.
/// Otherwise return zero.
fn transit(lon1: f64, lon2: f64) -> i32 {
    // Compute lon12 the same way as Geodesic::inverse.
    let (lon12, _) = ang_diff(lon1, lon2);
    let lon1 = ang_normalize(lon1);
    let lon2 = ang_normalize(lon2);
    if lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)) {
        1
    } else if lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0 {
        -1
    } else {
        0
    }
}

/// Compute exactly the parity of
///   `int(floor(lon2 / 360)) - int(floor(lon1 / 360))`
fn transitdirect(lon1: f64, lon2: f64) -> i32 {
    let lon1 = remainder(lon1, 2.0 * TD);
    let lon2 = remainder(lon2, 2.0 * TD);
    let outside = |lon: f64| if (0.0..TD).contains(&lon) { 0 } else { 1 };
    outside(lon2) - outside(lon1)
}

fn area_reduce_acc(
    mut area: Accumulator,
    area0: f64,
    crossings: i32,
    reverse: bool,
    sign: bool,
) -> f64 {
    area.remainder(area0);
    if crossings & 1 != 0 {
        area.add((if area.sum() < 0.0 { 1.0 } else { -1.0 }) * area0 / 2.0);
    }
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if !reverse {
        area.negate();
    }
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if sign {
        if area.sum() > area0 / 2.0 {
            area.add(-area0);
        } else if area.sum() <= -area0 / 2.0 {
            area.add(area0);
        }
    } else if area.sum() >= area0 {
        area.add(-area0);
    } else if area.sum() < 0.0 {
        area.add(area0);
    }
    0.0 + area.sum()
}

fn area_reduce(area: f64, area0: f64, crossings: i32, reverse: bool, sign: bool) -> f64 {
    let mut area = remainder(area, area0);
    if crossings & 1 != 0 {
        area += (if area < 0.0 { 1.0 } else { -1.0 }) * area0 / 2.0;
    }
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if !reverse {
        area *= -1.0;
    }
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if sign {
        if area > area0 / 2.0 {
            area -= area0;
        } else if area <= -area0 / 2.0 {
            area += area0;
        }
    } else if area >= area0 {
        area -= area0;
    } else if area < 0.0 {
        area += area0;
    }
    0.0 + area
}
