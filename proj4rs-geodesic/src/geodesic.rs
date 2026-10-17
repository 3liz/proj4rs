//! Ellipsoid and solution of the direct and inverse geodesic problems.
//!
//! This is a transcription of `geodesic.c` from GeographicLib/Proj.
//! Most of the internal comments have been retained, and the structure of
//! the computations is kept close to the original so that fixes in the
//! reference implementation can be easily ported.
// Comparisons are kept as in the reference implementation so that
// NaNs propagate the same way.
#![allow(clippy::neg_cmp_op_on_partial_ord)]
use crate::line::GeodesicLine;
use crate::math::*;
use std::f64::consts::PI;
use std::ops::{BitAnd, BitOr, BitOrAssign};

const GEODESIC_ORDER: usize = 6;
const NA1: usize = GEODESIC_ORDER;
pub(crate) const NC1: usize = GEODESIC_ORDER;
pub(crate) const NC1P: usize = GEODESIC_ORDER;
const NA2: usize = GEODESIC_ORDER;
pub(crate) const NC2: usize = GEODESIC_ORDER;
const NA3: usize = GEODESIC_ORDER;
const NA3X: usize = NA3;
pub(crate) const NC3: usize = GEODESIC_ORDER;
const NC3X: usize = (NC3 * (NC3 - 1)) / 2;
pub(crate) const NC4: usize = GEODESIC_ORDER;
const NC4X: usize = (NC4 * (NC4 + 1)) / 2;
pub(crate) const NC: usize = GEODESIC_ORDER + 1;

/// Bit mask specifying the quantities to compute (and, for a
/// [`GeodesicLine`], the capabilities of the line).
///
/// Masks may be combined with `|`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Mask(pub(crate) u32);

impl Mask {
    pub(crate) const CAP_C1: u32 = 1 << 0;
    pub(crate) const CAP_C1P: u32 = 1 << 1;
    pub(crate) const CAP_C2: u32 = 1 << 2;
    pub(crate) const CAP_C3: u32 = 1 << 3;
    pub(crate) const CAP_C4: u32 = 1 << 4;
    pub(crate) const OUT_ALL: u32 = 0x7F80;

    /// Calculate nothing
    pub const NONE: Self = Self(0);
    /// Calculate latitude
    pub const LATITUDE: Self = Self(1 << 7);
    /// Calculate longitude
    pub const LONGITUDE: Self = Self(1 << 8 | Self::CAP_C3);
    /// Calculate azimuth
    pub const AZIMUTH: Self = Self(1 << 9);
    /// Calculate distance
    pub const DISTANCE: Self = Self(1 << 10 | Self::CAP_C1);
    /// Allow distance as input
    pub const DISTANCE_IN: Self = Self(1 << 11 | Self::CAP_C1 | Self::CAP_C1P);
    /// Calculate reduced length
    pub const REDUCED_LENGTH: Self = Self(1 << 12 | Self::CAP_C1 | Self::CAP_C2);
    /// Calculate geodesic scale
    pub const GEODESIC_SCALE: Self = Self(1 << 13 | Self::CAP_C1 | Self::CAP_C2);
    /// Calculate area
    pub const AREA: Self = Self(1 << 14 | Self::CAP_C4);
    /// Calculate everything
    pub const ALL: Self = Self(Self::OUT_ALL | 0x1F);

    /// Return true if `self` contains all the output bits of `other`
    pub const fn contains(self, other: Self) -> bool {
        let out = other.0 & Self::OUT_ALL;
        self.0 & out == out
    }

    /// Return true if any of the bits of `bits` are set
    #[inline]
    pub(crate) const fn any(self, bits: u32) -> bool {
        self.0 & bits != 0
    }

    /// Return the output bits of `self`
    #[inline]
    pub(crate) const fn out(self) -> u32 {
        self.0 & Self::OUT_ALL
    }
}

impl BitOr for Mask {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for Mask {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0
    }
}

impl BitAnd for Mask {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self {
        Self(self.0 & rhs.0)
    }
}

/// Flags for the direct problem and the positions along a [`GeodesicLine`].
///
/// Flags may be combined with `|`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Flags(u32);

impl Flags {
    /// No flags
    pub const NONE: Self = Self(0);
    /// Position given in terms of arc distance (degrees) instead of meters
    pub const ARC_MODE: Self = Self(1 << 0);
    /// Unroll the longitude
    pub const LONG_UNROLL: Self = Self(1 << 15);

    /// Return true if all the bits of `other` are set
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl BitOr for Flags {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// Result of the general inverse problem.
///
/// Quantities not requested are set to NaN.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Inverse {
    /// Arc length from point 1 to point 2 \[degrees\]
    pub a12: f64,
    /// Distance from point 1 to point 2 \[meters\]
    pub s12: f64,
    /// Azimuth at point 1 \[degrees\]
    pub azi1: f64,
    /// (Forward) azimuth at point 2 \[degrees\]
    pub azi2: f64,
    /// Reduced length of the geodesic, `m12` \[meters\]
    pub m12: f64,
    /// Geodesic scale of point 2 relative to point 1, `M12` (dimensionless)
    pub scale12: f64,
    /// Geodesic scale of point 1 relative to point 2, `M21` (dimensionless)
    pub scale21: f64,
    /// Area under the geodesic, `S12` \[meters²\]
    pub area12: f64,
}

/// Result of the general direct problem.
///
/// Quantities not requested (or not supported by the
/// capabilities of a [`GeodesicLine`]) are set to NaN.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Position {
    /// Arc length from point 1 to point 2 \[degrees\]
    pub a12: f64,
    /// Latitude of point 2 \[degrees\]
    pub lat2: f64,
    /// Longitude of point 2 \[degrees\]
    pub lon2: f64,
    /// (Forward) azimuth at point 2 \[degrees\]
    pub azi2: f64,
    /// Distance from point 1 to point 2 \[meters\]
    pub s12: f64,
    /// Reduced length of the geodesic, `m12` \[meters\]
    pub m12: f64,
    /// Geodesic scale of point 2 relative to point 1, `M12` (dimensionless)
    pub scale12: f64,
    /// Geodesic scale of point 1 relative to point 2, `M21` (dimensionless)
    pub scale21: f64,
    /// Area under the geodesic, `S12` \[meters²\]
    pub area12: f64,
}

impl Default for Position {
    fn default() -> Self {
        Self {
            a12: f64::NAN,
            lat2: f64::NAN,
            lon2: f64::NAN,
            azi2: f64::NAN,
            s12: f64::NAN,
            m12: f64::NAN,
            scale12: f64::NAN,
            scale21: f64::NAN,
            area12: f64::NAN,
        }
    }
}

/// Output of `lengths`
#[derive(Default)]
struct Lengths {
    s12b: f64,
    m12b: f64,
    m0: f64,
    mm12: f64,
    mm21: f64,
}

/// Output of `lambda12`
struct Lambda12 {
    lam12: f64,
    salp2: f64,
    calp2: f64,
    sig12: f64,
    ssig1: f64,
    csig1: f64,
    ssig2: f64,
    csig2: f64,
    eps: f64,
    domg12: f64,
    dlam12: f64,
}

/// Output of `inverse_start`
struct InverseStart {
    sig12: f64,
    salp1: f64,
    calp1: f64,
    // Only valid if sig12 >= 0
    salp2: f64,
    calp2: f64,
    // Only valid for short lines
    dnm: f64,
}

/// Output of `gen_inverse_int`
pub(crate) struct InverseInt {
    pub a12: f64,
    pub s12: f64,
    pub salp1: f64,
    pub calp1: f64,
    pub salp2: f64,
    pub calp2: f64,
    pub m12: f64,
    pub mm12: f64,
    pub mm21: f64,
    pub ss12: f64,
}

/// Ellipsoid on which geodesic calculations are computed
#[derive(Clone)]
pub struct Geodesic {
    /// Semi-major axis
    pub(crate) a: f64,
    /// Flattening
    pub(crate) f: f64,
    pub(crate) f1: f64,
    pub(crate) e2: f64,
    pub(crate) ep2: f64,
    pub(crate) n: f64,
    pub(crate) b: f64,
    pub(crate) c2: f64,
    pub(crate) etol2: f64,
    a3x: [f64; NA3X],
    c3x: [f64; NC3X],
    c4x: [f64; NC4X],
}

impl std::fmt::Display for Geodesic {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Geodesic {{ a: {}, f: {} }}", self.a, self.f)
    }
}

impl std::fmt::Debug for Geodesic {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Geodesic {{ a: {}, f: {} }}", self.a, self.f)
    }
}

impl Geodesic {
    /// Create new Ellipsoid with semi-major axis `a` in meters and a flattening `f`
    ///
    /// ```rust
    /// use proj4rs_geodesic::Geodesic;
    /// let g = Geodesic::new(6_378_145.0, 1.0/298.25);
    /// println!("{}", g);
    /// // Geodesic { a: 6378145, f: 0.003352891869237217 }
    /// ```
    pub fn new(a: f64, f: f64) -> Self {
        let f1 = 1.0 - f;
        let e2 = f * (2.0 - f);
        let ep2 = e2 / sq(f1); // e2 / (1 - e2)
        let n = f / (2.0 - f);
        let b = a * f1;
        // authalic radius squared
        let c2 = (sq(a)
            + sq(b)
                * (if e2 == 0.0 {
                    1.0
                } else {
                    (if e2 > 0.0 {
                        atanh(e2.sqrt())
                    } else {
                        (-e2).sqrt().atan()
                    }) / e2.abs().sqrt()
                }))
            / 2.0;
        // The sig12 threshold for "really short".  Using the auxiliary sphere
        // solution with dnm computed at (bet1 + bet2) / 2, the relative error in
        // the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
        // (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a
        // given f and sig12, the max error occurs for lines near the pole.  If
        // the old rule for computing dnm = (dn1 + dn2)/2 is used, then the error
        // increases by a factor of 2.)  Setting this equal to epsilon gives
        // sig12 = etol2.  Here 0.1 is a safety factor (error decreased by 100)
        // and max(0.001, abs(f)) stops etol2 getting too large in the nearly
        // spherical case.
        let etol2 = 0.1 * TOL2 / (f.abs().max(0.001) * (1.0 - f / 2.0).min(1.0) / 2.0).sqrt();

        Self {
            a,
            f,
            f1,
            e2,
            ep2,
            n,
            b,
            c2,
            etol2,
            a3x: a3_coeff(n),
            c3x: c3_coeff(n),
            c4x: c4_coeff(n),
        }
    }

    /// Create the WGS84 ellipsoid
    pub fn wgs84() -> Self {
        const A: f64 = 6_378_137.0;
        const F: f64 = 1.0 / 298.257_223_563;
        Self::new(A, F)
    }

    /// The equatorial radius \[meters\]
    #[inline(always)]
    pub fn a(&self) -> f64 {
        self.a
    }

    /// The flattening
    #[inline(always)]
    pub fn f(&self) -> f64 {
        self.f
    }

    /// Compute a new location (`lat2`,`lon2`) from (`lat1`,`lon1`) a distance `s12` at an azimuth of `azi1`
    ///
    /// # Arguments
    ///   - lat1 - Latitude of 1st point [degrees] [-90.,90.]
    ///   - lon1 - Longitude of 1st point [degrees] [-180., 180.]
    ///   - azi1 - Azimuth at 1st point [degrees] [-180., 180.]
    ///   - s12 - Distance from 1st to 2nd point [meters] Value may be negative
    ///
    /// # Returns
    ///   - lat2 - Latitude of 2nd point [degrees]
    ///   - lon2 - Longitude of 2nd point [degrees]
    ///   - azi2 - Azimuth at 2nd point
    ///
    /// If either point is at a pole, the azimuth is defined by keeping the
    ///  longitude fixed, writing lat = ±(90° − ε), and taking the limit ε → 0+.
    ///  An arc length greater that 180° signifies a geodesic which is not a
    ///  shortest path. (For a prolate ellipsoid, an additional condition is
    ///  necessary for a shortest path: the longitudinal extent must not
    ///  exceed of 180°.)
    ///
    /// ```rust
    /// // Example, determine the point 10000 km NE of JFK:
    /// use proj4rs_geodesic::Geodesic;
    /// let g = Geodesic::wgs84();
    /// let (lat,lon,az) = g.direct(40.64, -73.78, 45.0, 10e6);
    /// assert_eq!(lat, 32.621100463725796);
    /// assert_eq!(lon, 49.05248709295982);
    /// assert_eq!(az,  140.40598587680074);
    /// ```
    ///
    pub fn direct(&self, lat1: f64, lon1: f64, azi1: f64, s12: f64) -> (f64, f64, f64) {
        let p = self.gen_direct(
            lat1,
            lon1,
            azi1,
            Flags::NONE,
            s12,
            Mask::LATITUDE | Mask::LONGITUDE | Mask::AZIMUTH,
        );
        (p.lat2, p.lon2, p.azi2)
    }

    /// The general direct geodesic problem.
    ///
    /// # Arguments
    ///   - lat1 - Latitude of 1st point [degrees] [-90.,90.]
    ///   - lon1 - Longitude of 1st point [degrees]
    ///   - azi1 - Azimuth at 1st point [degrees]
    ///   - flags - if `flags` contains [`Flags::ARC_MODE`], `s12_a12` is the
    ///     arc length (degrees) instead of the distance (meters);
    ///     [`Flags::LONG_UNROLL`] "unrolls" `lon2`, i.e. `lon2 - lon1` indicates
    ///     how many times and in what sense the geodesic encircles the ellipsoid.
    ///   - s12_a12 - Distance [meters] or arc length [degrees] from 1st to
    ///     2nd point; it can be negative.
    ///   - outmask - the quantities to compute.
    ///
    /// The values of `lon2` and `azi2` returned are in the range [-180°, 180°].
    pub fn gen_direct(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        flags: Flags,
        s12_a12: f64,
        outmask: Mask,
    ) -> Position {
        // Automatically supply DISTANCE_IN if necessary
        let caps = if flags.contains(Flags::ARC_MODE) {
            outmask
        } else {
            outmask | Mask::DISTANCE_IN
        };
        GeodesicLine::new(self, lat1, lon1, azi1, caps).gen_position(flags, s12_a12, outmask)
    }

    /// Compute distance and azimuth from (`lat1`,`lon1`) to (`lat2`,`lon2`)
    ///
    /// # Arguments
    ///   - lat1: Latitude of 1st point [degrees] [-90., 90.]
    ///   - lon1: Longitude of 1st point [degrees] [-180., 180.]
    ///   - lat2: Latitude of 2nd point [degrees] [-90. 90]
    ///   - lon2: Longitude of 2nd point [degrees] [-180., 180.]
    ///
    /// # Returns
    ///   - s12: Distance from 1st to 2nd point [meters]
    ///   - azi1: Azimuth at 1st point [degrees]
    ///   - azi2: Azimuth at 2nd point [degrees]
    ///
    /// If either point is at a pole, the azimuth is defined by keeping the
    ///   longitude fixed, writing lat = ±(90° − ε), and taking the limit ε → 0+.
    ///
    /// The solution to the inverse problem is found using Newton's method.
    ///  If this fails to converge (this is very unlikely in geodetic applications
    ///  but does occur for very eccentric ellipsoids), then the bisection method
    ///  is used to refine the solution.
    ///
    /// ```rust
    /// // Example, determine the distance between JFK and Singapore Changi Airport:
    /// use proj4rs_geodesic::Geodesic;
    /// let g = Geodesic::wgs84();
    /// let (jfk_lat, jfk_lon) = (40.64, -73.78);
    /// let (sin_lat, sin_lon) = (1.36, 103.99);
    /// let (m, a1, a2) = g.inverse(jfk_lat, jfk_lon, sin_lat, sin_lon);
    /// assert_eq!(m,  15347512.94051294);  // Distance meters
    /// assert_eq!(a1, 3.3057734780176125); // Azimuth at 1st point
    /// assert_eq!(a2, 177.48784020815515); // Azimuth at 2nd point (forward)
    /// ```
    ///
    pub fn inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> (f64, f64, f64) {
        let r = self.gen_inverse(lat1, lon1, lat2, lon2, Mask::DISTANCE);
        (r.s12, r.azi1, r.azi2)
    }

    /// The general inverse geodesic calculation.
    ///
    /// # Arguments
    ///   - lat1: Latitude of 1st point [degrees] [-90., 90.]
    ///   - lon1: Longitude of 1st point [degrees]
    ///   - lat2: Latitude of 2nd point [degrees] [-90. 90]
    ///   - lon2: Longitude of 2nd point [degrees]
    ///   - outmask: the quantities to compute, the arc length and the
    ///     azimuths are always computed.
    pub fn gen_inverse(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        outmask: Mask,
    ) -> Inverse {
        let r = self.gen_inverse_int(lat1, lon1, lat2, lon2, outmask);
        let get = |m: Mask, v: f64| if outmask.any(m.out()) { v } else { f64::NAN };
        Inverse {
            a12: r.a12,
            s12: get(Mask::DISTANCE, r.s12),
            azi1: atan2d(r.salp1, r.calp1),
            azi2: atan2d(r.salp2, r.calp2),
            m12: get(Mask::REDUCED_LENGTH, r.m12),
            scale12: get(Mask::GEODESIC_SCALE, r.mm12),
            scale21: get(Mask::GEODESIC_SCALE, r.mm21),
            area12: get(Mask::AREA, r.ss12),
        }
    }

    /// Create a [`GeodesicLine`] starting at (`lat1`, `lon1`) with azimuth `azi1`.
    ///
    /// `caps` specifies the capabilities of the line, see [`GeodesicLine::new`].
    pub fn line(&self, lat1: f64, lon1: f64, azi1: f64, caps: Mask) -> GeodesicLine {
        GeodesicLine::new(self, lat1, lon1, azi1, caps)
    }

    /// Create a [`GeodesicLine`] in terms of the direct geodesic problem.
    ///
    /// Point 3 of the line is set to point 2 of the direct problem.
    pub fn direct_line(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        s12: f64,
        caps: Mask,
    ) -> GeodesicLine {
        self.gen_direct_line(lat1, lon1, azi1, Flags::NONE, s12, caps)
    }

    /// Create a [`GeodesicLine`] in terms of the direct geodesic problem
    /// specified in terms of either distance (meters) or arc length (degrees)
    /// depending on `flags`.
    ///
    /// Point 3 of the line is set to point 2 of the direct problem.
    pub fn gen_direct_line(
        &self,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        flags: Flags,
        s12_a12: f64,
        caps: Mask,
    ) -> GeodesicLine {
        let mut l = GeodesicLine::new(self, lat1, lon1, azi1, caps);
        l.gen_set_distance(flags, s12_a12);
        l
    }

    /// Create a [`GeodesicLine`] in terms of the inverse geodesic problem.
    ///
    /// Point 3 of the line is set to point 2 of the inverse problem.
    ///
    /// ```rust
    /// // Compute way points between JFK and Singapore Changi Airport
    /// use proj4rs_geodesic::{Geodesic, Mask};
    /// let g = Geodesic::wgs84();
    /// let l = g.inverse_line(40.64, -73.78, 1.36, 103.99, Mask::NONE);
    /// let waypoints: Vec<_> = (0..=100)
    ///     .map(|i| l.position(i as f64 * l.s13() * 0.01))
    ///     .collect();
    /// let (lat, lon, _) = waypoints[100];
    /// assert!((lat - 1.36).abs() < 1e-12 && (lon - 103.99).abs() < 1e-12);
    /// ```
    pub fn inverse_line(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        caps: Mask,
    ) -> GeodesicLine {
        let r = self.gen_inverse_int(lat1, lon1, lat2, lon2, Mask::NONE);
        let azi1 = atan2d(r.salp1, r.calp1);
        let mut caps = if caps.0 != 0 {
            caps
        } else {
            Mask::DISTANCE_IN | Mask::LONGITUDE
        };
        // Ensure that a12 can be converted to a distance
        if caps.any(Mask::DISTANCE_IN.out()) {
            caps |= Mask::DISTANCE;
        }
        let mut l = GeodesicLine::new_int(self, lat1, lon1, azi1, r.salp1, r.calp1, caps);
        l.set_arc(r.a12);
        l
    }

    #[allow(clippy::approx_constant)]
    pub(crate) fn gen_inverse_int(
        &self,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        outmask: Mask,
    ) -> InverseInt {
        let outmask = outmask & Mask(Mask::OUT_ALL);
        let want_scale = outmask.any(Mask::GEODESIC_SCALE.out());

        let (mut s12x, mut m12x) = (0.0, 0.0);
        let (mut mm12, mut mm21, mut ss12) = (0.0, 0.0, 0.0);
        let mut a12 = 0.0;
        let mut sig12;
        let (mut salp1, mut calp1, mut salp2, mut calp2) = (0.0, 0.0, 0.0, 0.0);
        // somg12 == 2 marks that it needs to be calculated
        let (mut omg12, mut somg12, mut comg12) = (0.0, 2.0, 0.0);

        // Compute longitude difference (ang_diff does this carefully).  Result is
        // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
        // east-going and meridional geodesics.
        let (mut lon12, mut lon12s) = ang_diff(lon1, lon2);
        // Make longitude difference positive.
        let mut lonsign: f64 = if lon12.is_sign_negative() { -1.0 } else { 1.0 };
        lon12 *= lonsign;
        lon12s *= lonsign;
        let lam12 = lon12 * DEGREE;
        // Calculate sincos of lon12 + error (this applies ang_round internally).
        let (slam12, clam12) = sincosde(lon12, lon12s);
        // the supplementary longitude difference
        let lon12s = (HD - lon12) - lon12s;

        // If really close to the equator, treat as on equator.
        let mut lat1 = ang_round(lat_fix(lat1));
        let mut lat2 = ang_round(lat_fix(lat2));
        // Swap points so that point with higher (abs) latitude is point 1
        // If one latitude is a nan, then it becomes lat1.
        let swapp: f64 = if lat1.abs() < lat2.abs() || lat2.is_nan() {
            -1.0
        } else {
            1.0
        };
        if swapp < 0.0 {
            lonsign *= -1.0;
            std::mem::swap(&mut lat1, &mut lat2);
        }
        // Make lat1 <= -0
        let latsign: f64 = if lat1.is_sign_negative() { 1.0 } else { -1.0 };
        lat1 *= latsign;
        lat2 *= latsign;
        // Now we have
        //
        //     0 <= lon12 <= 180
        //     -90 <= lat1 <= -0
        //     lat1 <= lat2 <= -lat1
        //
        // longsign, swapp, latsign register the transformation to bring the
        // coordinates to this canonical form.  In all cases, 1 means no change
        // was made.  We make these transformations so that there are few cases
        // to check, e.g., on verifying quadrants in atan2.  In addition, this
        // enforces some symmetries in the results returned.

        let (sbet1, cbet1) = sincosd(lat1);
        // Ensure cbet1 = +epsilon at poles
        let (sbet1, cbet1) = norm2(sbet1 * self.f1, cbet1);
        let cbet1 = cbet1.max(TINY);

        let (sbet2, cbet2) = sincosd(lat2);
        // Ensure cbet2 = +epsilon at poles
        let (mut sbet2, mut cbet2) = norm2(sbet2 * self.f1, cbet2);
        cbet2 = cbet2.max(TINY);

        // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
        // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1
        // is a better measure.  This logic is used in assigning calp2 in
        // lambda12.  Sometimes these quantities vanish and in that case we force
        // bet2 = +/- bet1 exactly.  An example where is is necessary is the
        // inverse problem
        // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
        // which failed with Visual Studio 10 (Release and Debug)
        if cbet1 < -sbet1 {
            if cbet2 == cbet1 {
                sbet2 = sbet1.copysign(sbet2);
            }
        } else if sbet2.abs() == -sbet1 {
            cbet2 = cbet1;
        }

        let dn1 = (1.0 + self.ep2 * sq(sbet1)).sqrt();
        let dn2 = (1.0 + self.ep2 * sq(sbet2)).sqrt();

        let lengths_mask = if want_scale {
            Mask::DISTANCE | Mask::REDUCED_LENGTH | Mask::GEODESIC_SCALE
        } else {
            Mask::DISTANCE | Mask::REDUCED_LENGTH
        };

        let mut meridian = lat1 == -QD || slam12 == 0.0;

        if meridian {
            // Endpoints are on a single full meridian, so the geodesic might lie
            // on a meridian.

            calp1 = clam12;
            salp1 = slam12; // Head to the target longitude
            calp2 = 1.0;
            salp2 = 0.0; // At the target we're heading north

            // tan(bet) = tan(sig) * cos(alp)
            let (ssig1, csig1) = (sbet1, calp1 * cbet1);
            let (ssig2, csig2) = (sbet2, calp2 * cbet2);

            // sig12 = sig2 - sig1
            sig12 = ((csig1 * ssig2 - ssig1 * csig2).max(0.0) + 0.0)
                .atan2(csig1 * csig2 + ssig1 * ssig2);
            let l = self.lengths(
                self.n,
                sig12,
                ssig1,
                csig1,
                dn1,
                ssig2,
                csig2,
                dn2,
                cbet1,
                cbet2,
                lengths_mask,
            );
            s12x = l.s12b;
            m12x = l.m12b;
            if want_scale {
                mm12 = l.mm12;
                mm21 = l.mm21;
            }
            // Add the check for sig12 since zero length geodesics might yield
            // m12 < 0.  Test case was
            //
            //    echo 20.001 0 20.001 0 | GeodSolve -i
            //
            // In fact, we will have sig12 > pi/2 for meridional geodesic which
            // is not a shortest path.
            if sig12 < 1.0 || m12x >= 0.0 {
                // Need at least 2, to handle 90 0 90 180
                if sig12 < 3.0 * TINY ||
                    // Prevent negative s12 or m12 for short lines
                    (sig12 < TOL0 && (s12x < 0.0 || m12x < 0.0))
                {
                    sig12 = 0.0;
                    m12x = 0.0;
                    s12x = 0.0;
                }
                m12x *= self.b;
                s12x *= self.b;
                a12 = sig12 / DEGREE;
            } else {
                // m12 < 0, i.e., prolate and too close to anti-podal
                meridian = false;
            }
        }

        if !meridian && sbet1 == 0.0 && // and sbet2 == 0
            // Mimic the way lambda12 works with calp1 = 0
            (self.f <= 0.0 || lon12s >= self.f * HD)
        {
            // Geodesic runs along equator
            calp1 = 0.0;
            calp2 = 0.0;
            salp1 = 1.0;
            salp2 = 1.0;
            s12x = self.a * lam12;
            sig12 = lam12 / self.f1;
            omg12 = sig12;
            m12x = self.b * sig12.sin();
            if want_scale {
                mm12 = sig12.cos();
                mm21 = mm12;
            }
            a12 = lon12 / self.f1;
        } else if !meridian {
            // Now point1 and point2 belong within a hemisphere bounded by a
            // meridian and geodesic is neither meridional or equatorial.

            // Figure a starting point for Newton's method
            let start =
                self.inverse_start(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12);
            sig12 = start.sig12;
            salp1 = start.salp1;
            calp1 = start.calp1;

            if sig12 >= 0.0 {
                // Short lines (inverse_start sets salp2, calp2, dnm)
                let dnm = start.dnm;
                salp2 = start.salp2;
                calp2 = start.calp2;
                s12x = sig12 * self.b * dnm;
                m12x = sq(dnm) * self.b * (sig12 / dnm).sin();
                if want_scale {
                    mm12 = (sig12 / dnm).cos();
                    mm21 = mm12;
                }
                a12 = sig12 / DEGREE;
                omg12 = lam12 / (self.f1 * dnm);
            } else {
                // Newton's method.  This is a straightforward solution of f(alp1)
                // = lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has
                // exactly one root in the interval (0, pi) and its derivative is
                // positive at the root.  Thus f(alp) is positive for alp > alp1
                // and negative for alp < alp1.  During the course of the
                // iteration, a range (alp1a, alp1b) is maintained which brackets
                // the root and with each evaluation of f(alp) the range is
                // shrunk, if possible.  Newton's method is restarted whenever the
                // derivative of f is negative (because the new value of alp1 is
                // then further from the solution) or if the new estimate of alp1
                // lies outside (0,pi); in this case, the new starting guess is
                // taken to be (alp1a + alp1b) / 2.
                let (mut ssig1, mut csig1, mut ssig2, mut csig2, mut eps, mut domg12);
                // Bracketing range
                let (mut salp1a, mut calp1a, mut salp1b, mut calp1b) = (TINY, 1.0, TINY, -1.0);
                let mut tripn = false;
                let mut tripb = false;
                let mut numit = 0;
                loop {
                    // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
                    // WGS84 and random input: mean = 2.85, sd = 0.60
                    let r = self.lambda12(
                        sbet1,
                        cbet1,
                        dn1,
                        sbet2,
                        cbet2,
                        dn2,
                        salp1,
                        calp1,
                        slam12,
                        clam12,
                        numit < MAXIT1,
                    );
                    let v = r.lam12;
                    let dv = r.dlam12;
                    salp2 = r.salp2;
                    calp2 = r.calp2;
                    sig12 = r.sig12;
                    ssig1 = r.ssig1;
                    csig1 = r.csig1;
                    ssig2 = r.ssig2;
                    csig2 = r.csig2;
                    eps = r.eps;
                    domg12 = r.domg12;

                    if tripb ||
                        // Reversed test to allow escape with NaNs
                        !(v.abs() >= (if tripn { 8.0 } else { 1.0 }) * TOL0) ||
                        // Enough bisections to get accurate result
                        numit == MAXIT2
                    {
                        break;
                    }
                    // Update bracketing values
                    if v > 0.0 && (numit > MAXIT1 || calp1 / salp1 > calp1b / salp1b) {
                        salp1b = salp1;
                        calp1b = calp1;
                    } else if v < 0.0 && (numit > MAXIT1 || calp1 / salp1 < calp1a / salp1a) {
                        salp1a = salp1;
                        calp1a = calp1;
                    }
                    let newton = numit < MAXIT1 && dv > 0.0;
                    numit += 1;
                    if newton {
                        let dalp1 = -v / dv;
                        if dalp1.abs() < PI {
                            let (sdalp1, cdalp1) = dalp1.sin_cos();
                            let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                            if nsalp1 > 0.0 {
                                (salp1, calp1) = norm2(nsalp1, calp1 * cdalp1 - salp1 * sdalp1);
                                // In some regimes we don't get quadratic
                                // convergence because slope -> 0.  So use
                                // convergence conditions based on epsilon
                                // instead of sqrt(epsilon).
                                tripn = v.abs() <= 16.0 * TOL0;
                                continue;
                            }
                        }
                    }
                    // Either dv was not positive or updated value was outside
                    // legal range.  Use the midpoint of the bracket as the next
                    // estimate.  This mechanism is not needed for the WGS84
                    // ellipsoid, but it does catch problems with more eccentric
                    // ellipsoids.  Its efficacy is such for the WGS84 test set
                    // with the starting guess set to alp1 = 90deg:
                    // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
                    // WGS84 and random input: mean = 4.74, sd = 0.99
                    (salp1, calp1) = norm2((salp1a + salp1b) / 2.0, (calp1a + calp1b) / 2.0);
                    tripn = false;
                    tripb = (salp1a - salp1).abs() + (calp1a - calp1) < TOLB
                        || (salp1 - salp1b).abs() + (calp1 - calp1b) < TOLB;
                }
                let l = self.lengths(
                    eps,
                    sig12,
                    ssig1,
                    csig1,
                    dn1,
                    ssig2,
                    csig2,
                    dn2,
                    cbet1,
                    cbet2,
                    lengths_mask,
                );
                s12x = l.s12b * self.b;
                m12x = l.m12b * self.b;
                if want_scale {
                    mm12 = l.mm12;
                    mm21 = l.mm21;
                }
                a12 = sig12 / DEGREE;
                if outmask.any(Mask::AREA.out()) {
                    // omg12 = lam12 - domg12
                    let (sdomg12, cdomg12) = domg12.sin_cos();
                    somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                    comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                }
            }
        }

        // Convert -0 to 0
        let s12 = if outmask.any(Mask::DISTANCE.out()) {
            0.0 + s12x
        } else {
            0.0
        };
        let m12 = if outmask.any(Mask::REDUCED_LENGTH.out()) {
            0.0 + m12x
        } else {
            0.0
        };

        if outmask.any(Mask::AREA.out()) {
            // From lambda12: sin(alp1) * cos(bet1) = sin(alp0)
            let salp0 = salp1 * cbet1;
            let calp0 = calp1.hypot(salp1 * sbet1); // calp0 > 0
            if calp0 != 0.0 && salp0 != 0.0 {
                // From lambda12: tan(bet) = tan(sig) * cos(alp)
                let (ssig1, csig1) = norm2(sbet1, calp1 * cbet1);
                let (ssig2, csig2) = norm2(sbet2, calp2 * cbet2);
                let k2 = sq(calp0) * self.ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                let a4 = sq(self.a) * calp0 * salp0 * self.e2;
                let c4a = self.c4f(eps);
                let b41 = cos_series(ssig1, csig1, &c4a);
                let b42 = cos_series(ssig2, csig2, &c4a);
                ss12 = a4 * (b42 - b41);
            } else {
                // Avoid problems with indeterminate sig1, sig2 on equator
                ss12 = 0.0;
            }

            if !meridian && somg12 == 2.0 {
                (somg12, comg12) = omg12.sin_cos();
            }

            let alp12 = if !meridian &&
                // omg12 < 3/4 * pi
                comg12 > -0.7071 &&     // Long difference not too big
                sbet2 - sbet1 < 1.75
            {
                // Lat difference not too big
                // Use tan(Gamma/2) = tan(omg12/2)
                // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
                // with tan(x/2) = sin(x)/(1+cos(x))
                let (domg12, dbet1, dbet2) = (1.0 + comg12, 1.0 + cbet1, 1.0 + cbet2);
                2.0 * (somg12 * (sbet1 * dbet2 + sbet2 * dbet1))
                    .atan2(domg12 * (sbet1 * sbet2 + dbet1 * dbet2))
            } else {
                // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                let mut salp12 = salp2 * calp1 - calp2 * salp1;
                let mut calp12 = calp2 * calp1 + salp2 * salp1;
                // The right thing appears to happen if alp1 = +/-180 and alp2 =
                // 0, viz salp12 = -0 and alp12 = -180.  However this depends on
                // the sign being attached to 0 correctly.  The following ensures
                // the correct behavior.
                if salp12 == 0.0 && calp12 < 0.0 {
                    salp12 = TINY * calp1;
                    calp12 = -1.0;
                }
                salp12.atan2(calp12)
            };
            ss12 += self.c2 * alp12;
            ss12 *= swapp * lonsign * latsign;
            // Convert -0 to 0
            ss12 += 0.0;
        }

        // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
        if swapp < 0.0 {
            std::mem::swap(&mut salp1, &mut salp2);
            std::mem::swap(&mut calp1, &mut calp2);
            if want_scale {
                std::mem::swap(&mut mm12, &mut mm21);
            }
        }

        salp1 *= swapp * lonsign;
        calp1 *= swapp * latsign;
        salp2 *= swapp * lonsign;
        calp2 *= swapp * latsign;

        // Returned value in [0, 180]
        InverseInt {
            a12,
            s12,
            salp1,
            calp1,
            salp2,
            calp2,
            m12,
            mm12,
            mm21,
            ss12,
        }
    }

    /// Return `m12b = (reduced length)/b`; also calculate `s12b = distance/b`,
    /// and `m0` = coefficient of secular term in expression for reduced length.
    ///
    /// `s12b` is computed if `outmask` contains [`Mask::DISTANCE`], `m12b` and
    /// `m0` if it contains [`Mask::REDUCED_LENGTH`] and `M12`, `M21` if it
    /// contains [`Mask::GEODESIC_SCALE`].
    #[allow(clippy::too_many_arguments)]
    fn lengths(
        &self,
        eps: f64,
        sig12: f64,
        ssig1: f64,
        csig1: f64,
        dn1: f64,
        ssig2: f64,
        csig2: f64,
        dn2: f64,
        cbet1: f64,
        cbet2: f64,
        outmask: Mask,
    ) -> Lengths {
        let distance = outmask.any(Mask::DISTANCE.out());
        let reduced = outmask.any(Mask::REDUCED_LENGTH.out());
        let scale = outmask.any(Mask::GEODESIC_SCALE.out());

        let mut r = Lengths::default();
        let (mut a1, mut a2, mut j12) = (0.0, 0.0, 0.0);
        let mut ca = [0.0; NC];
        let mut cb = [0.0; NC];

        let redlp = reduced || scale;
        if distance || redlp {
            a1 = a1m1f(eps);
            ca = c1f(eps);
            if redlp {
                a2 = a2m1f(eps);
                cb = c2f(eps);
                r.m0 = a1 - a2;
                a2 += 1.0;
            }
            a1 += 1.0;
        }
        if distance {
            let b1 =
                sin_series(ssig2, csig2, &ca[1..=NC1]) - sin_series(ssig1, csig1, &ca[1..=NC1]);
            // Missing a factor of b
            r.s12b = a1 * (sig12 + b1);
            if redlp {
                let b2 =
                    sin_series(ssig2, csig2, &cb[1..=NC2]) - sin_series(ssig1, csig1, &cb[1..=NC2]);
                j12 = r.m0 * sig12 + (a1 * b1 - a2 * b2);
            }
        } else if redlp {
            // Assume here that NC1 >= NC2
            for l in 1..=NC2 {
                cb[l] = a1 * ca[l] - a2 * cb[l];
            }
            j12 = r.m0 * sig12
                + (sin_series(ssig2, csig2, &cb[1..=NC2]) - sin_series(ssig1, csig1, &cb[1..=NC2]));
        }
        if reduced {
            // Missing a factor of b.
            // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
            // accurate cancellation in the case of coincident points.
            r.m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * j12;
        }
        if scale {
            let csig12 = csig1 * csig2 + ssig1 * ssig2;
            let t = self.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
            r.mm12 = csig12 + (t * ssig2 - csig2 * j12) * ssig1 / dn1;
            r.mm21 = csig12 - (t * ssig1 - csig1 * j12) * ssig2 / dn2;
        }
        r
    }

    /// Return a starting point for Newton's method in salp1 and calp1
    /// (function value is -1).  If Newton's method doesn't need to be used,
    /// return also salp2 and calp2 and function value is sig12.
    #[allow(clippy::too_many_arguments)]
    fn inverse_start(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        lam12: f64,
        slam12: f64,
        clam12: f64,
    ) -> InverseStart {
        let (mut salp2, mut calp2, mut dnm) = (0.0, 0.0, 0.0);

        let mut sig12 = -1.0; // Return value
        // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
        let sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
        let cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
        let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
        let sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
        let (mut somg12, mut comg12) = if shortline {
            // sin((bet1+bet2)/2)^2
            // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
            let mut sbetm2 = sq(sbet1 + sbet2);
            sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
            dnm = (1.0 + self.ep2 * sbetm2).sqrt();
            let omg12 = lam12 / (self.f1 * dnm);
            omg12.sin_cos()
        } else {
            (slam12, clam12)
        };

        let mut salp1 = cbet2 * somg12;
        let mut calp1 = if comg12 >= 0.0 {
            sbet12 + cbet2 * sbet1 * sq(somg12) / (1.0 + comg12)
        } else {
            sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12)
        };

        let ssig12 = salp1.hypot(calp1);
        let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

        if shortline && ssig12 < self.etol2 {
            // really short lines
            salp2 = cbet1 * somg12;
            calp2 = sbet12
                - cbet1
                    * sbet2
                    * (if comg12 >= 0.0 {
                        sq(somg12) / (1.0 + comg12)
                    } else {
                        1.0 - comg12
                    });
            (salp2, calp2) = norm2(salp2, calp2);
            // Set return value
            sig12 = ssig12.atan2(csig12);
        } else if self.n.abs() > 0.1 || // No astroid calc if too eccentric
            csig12 >= 0.0 ||
            ssig12 >= 6.0 * self.n.abs() * PI * sq(cbet1)
        {
            // Nothing to do, zeroth order spherical approximation is OK
        } else {
            // Scale lam12 and bet2 to x, y coordinate system where antipodal
            // point is at origin and singular point is at y = 0, x = -1.
            let (x, y, lamscale);
            let lam12x = (-slam12).atan2(-clam12); // lam12 - pi
            if self.f >= 0.0 {
                // In fact f == 0 does not get here
                // x = dlong, y = dlat
                let k2 = sq(sbet1) * self.ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                lamscale = self.f * cbet1 * self.a3f(eps) * PI;
                let betscale = lamscale * cbet1;

                x = lam12x / lamscale;
                y = sbet12a / betscale;
            } else {
                // f < 0
                // x = dlat, y = dlong
                let cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
                let bet12a = sbet12a.atan2(cbet12a);
                // In the case of lon12 = 180, this repeats a calculation made
                // in Inverse.
                let l = self.lengths(
                    self.n,
                    PI + bet12a,
                    sbet1,
                    -cbet1,
                    dn1,
                    sbet2,
                    cbet2,
                    dn2,
                    cbet1,
                    cbet2,
                    Mask::REDUCED_LENGTH,
                );
                x = -1.0 + l.m12b / (cbet1 * cbet2 * l.m0 * PI);
                let betscale = if x < -0.01 {
                    sbet12a / x
                } else {
                    -self.f * sq(cbet1) * PI
                };
                lamscale = betscale / cbet1;
                y = lam12x / lamscale;
            }

            if y > -TOL1 && x > -1.0 - XTHRESH {
                // strip near cut
                if self.f >= 0.0 {
                    salp1 = (-x).min(1.0);
                    calp1 = -(1.0 - sq(salp1)).sqrt();
                } else {
                    calp1 = x.max(if x > -TOL1 { 0.0 } else { -1.0 });
                    salp1 = (1.0 - sq(calp1)).sqrt();
                }
            } else {
                // Estimate alp1, by solving the astroid problem.
                //
                // Could estimate alpha1 = theta + pi/2, directly, i.e.,
                //   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
                //   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
                //
                // However, it's better to estimate omg12 from astroid and use
                // spherical formula to compute alp1.  This reduces the mean
                // number of Newton iterations for astroid cases from 2.24 (min
                // 0, max 6) to 2.12 (min 0 max 5).  The changes in the number
                // of iterations are as follows:
                //
                // change percent
                //    1       5
                //    0      78
                //   -1      16
                //   -2       0.6
                //   -3       0.04
                //   -4       0.002
                //
                // The histogram of iterations is (m = number of iterations
                // estimating alp1 directly, n = number of iterations estimating
                // via omg12, total number of trials = 148605):
                //
                //  iter    m      n
                //    0   148    186
                //    1 13046  13845
                //    2 93315 102225
                //    3 36189  32341
                //    4  5396      7
                //    5   455      1
                //    6    56      0
                //
                // Because omg12 is near pi, estimate work with omg12a = pi - omg12
                let k = astroid(x, y);
                let omg12a = lamscale
                    * if self.f >= 0.0 {
                        -x * k / (1.0 + k)
                    } else {
                        -y * (1.0 + k) / k
                    };
                somg12 = omg12a.sin();
                comg12 = -omg12a.cos();
                // Update spherical estimate of alp1 using omg12 instead of lam12
                salp1 = cbet2 * somg12;
                calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12);
            }
        }
        // Sanity check on starting guess.  Backwards check allows NaN through.
        if !(salp1 <= 0.0) {
            (salp1, calp1) = norm2(salp1, calp1);
        } else {
            salp1 = 1.0;
            calp1 = 0.0;
        }

        InverseStart {
            sig12,
            salp1,
            calp1,
            salp2,
            calp2,
            dnm,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn lambda12(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        salp1: f64,
        mut calp1: f64,
        slam120: f64,
        clam120: f64,
        diffp: bool,
    ) -> Lambda12 {
        if sbet1 == 0.0 && calp1 == 0.0 {
            // Break degeneracy of equatorial line.  This case has already been
            // handled.
            calp1 = -TINY;
        }

        // sin(alp1) * cos(bet1) = sin(alp0)
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1); // calp0 > 0

        // tan(bet1) = tan(sig1) * cos(alp1)
        // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
        let somg1 = salp0 * sbet1;
        let comg1 = calp1 * cbet1;
        let (ssig1, csig1) = norm2(sbet1, comg1);
        // norm2(somg1, comg1); -- don't need to normalize!

        // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
        // about this case, since this can yield singularities in the Newton
        // iteration.
        // sin(alp2) * cos(bet2) = sin(alp0)
        let salp2 = if cbet2 != cbet1 { salp0 / cbet2 } else { salp1 };
        // calp2 = sqrt(1 - sq(salp2))
        //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
        // and subst for calp0 and rearrange to give (choose positive sqrt
        // to give alp2 in [0, pi/2]).
        let calp2 = if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
            (sq(calp1 * cbet1)
                + if cbet1 < -sbet1 {
                    (cbet2 - cbet1) * (cbet1 + cbet2)
                } else {
                    (sbet1 - sbet2) * (sbet1 + sbet2)
                })
            .sqrt()
                / cbet2
        } else {
            calp1.abs()
        };
        // tan(bet2) = tan(sig2) * cos(alp2)
        // tan(omg2) = sin(alp0) * tan(sig2).
        let somg2 = salp0 * sbet2;
        let comg2 = calp2 * cbet2;
        let (ssig2, csig2) = norm2(sbet2, comg2);
        // norm2(somg2, comg2); -- don't need to normalize!

        // sig12 = sig2 - sig1, limit to [0, pi]
        let sig12 =
            ((csig1 * ssig2 - ssig1 * csig2).max(0.0) + 0.0).atan2(csig1 * csig2 + ssig1 * ssig2);

        // omg12 = omg2 - omg1, limit to [0, pi]
        let somg12 = (comg1 * somg2 - somg1 * comg2).max(0.0) + 0.0;
        let comg12 = comg1 * comg2 + somg1 * somg2;
        // eta = omg12 - lam120
        let eta = (somg12 * clam120 - comg12 * slam120).atan2(comg12 * clam120 + somg12 * slam120);
        let k2 = sq(calp0) * self.ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
        let c3a = self.c3f(eps);
        let b312 = sin_series(ssig2, csig2, &c3a[1..NC3]) - sin_series(ssig1, csig1, &c3a[1..NC3]);
        let domg12 = -self.f * self.a3f(eps) * salp0 * (sig12 + b312);
        let lam12 = eta + domg12;

        let dlam12 = if !diffp {
            0.0
        } else if calp2 == 0.0 {
            -2.0 * self.f1 * dn1 / sbet1
        } else {
            let l = self.lengths(
                eps,
                sig12,
                ssig1,
                csig1,
                dn1,
                ssig2,
                csig2,
                dn2,
                cbet1,
                cbet2,
                Mask::REDUCED_LENGTH,
            );
            l.m12b * (self.f1 / (calp2 * cbet2))
        };

        Lambda12 {
            lam12,
            salp2,
            calp2,
            sig12,
            ssig1,
            csig1,
            ssig2,
            csig2,
            eps,
            domg12,
            dlam12,
        }
    }

    /// Evaluate A3
    pub(crate) fn a3f(&self, eps: f64) -> f64 {
        polyval(&self.a3x, eps)
    }

    /// Evaluate C3 coeffs
    ///
    /// Elements c[1] through c[NC3 - 1] are set
    pub(crate) fn c3f(&self, eps: f64) -> [f64; NC3] {
        let mut c = [0.0; NC3];
        let mut mult = 1.0;
        let mut o = 0;
        for (l, cl) in c.iter_mut().enumerate().skip(1) {
            let m = NC3 - l - 1; // order of polynomial in eps
            mult *= eps;
            *cl = mult * polyval(&self.c3x[o..=o + m], eps);
            o += m + 1;
        }
        c
    }

    /// Evaluate C4 coeffs
    ///
    /// Elements c[0] through c[NC4 - 1] are set
    pub(crate) fn c4f(&self, eps: f64) -> [f64; NC4] {
        let mut c = [0.0; NC4];
        let mut mult = 1.0;
        let mut o = 0;
        for (l, cl) in c.iter_mut().enumerate() {
            let m = NC4 - l - 1; // order of polynomial in eps
            *cl = mult * polyval(&self.c4x[o..=o + m], eps);
            o += m + 1;
            mult *= eps;
        }
        c
    }
}

/// Evaluate `sum(c[i] * sin(2 * (i+1) * x), i, 0, n-1)` using Clenshaw
/// summation, where `n = c.len()`.
pub(crate) fn sin_series(sinx: f64, cosx: f64, c: &[f64]) -> f64 {
    let (y0, _) = clenshaw(sinx, cosx, c);
    2.0 * sinx * cosx * y0 // sin(2 * x) * y0
}

/// Evaluate `sum(c[i] * cos((2*i+1) * x), i, 0, n-1)` using Clenshaw
/// summation, where `n = c.len()`.
pub(crate) fn cos_series(sinx: f64, cosx: f64, c: &[f64]) -> f64 {
    let (y0, y1) = clenshaw(sinx, cosx, c);
    cosx * (y0 - y1) // cos(x) * (y0 - y1)
}

/// Approx operation count = (n + 5) mult and (2 * n + 2) add
#[inline]
fn clenshaw(sinx: f64, cosx: f64, c: &[f64]) -> (f64, f64) {
    let ar = 2.0 * (cosx - sinx) * (cosx + sinx); // 2 * cos(2 * x)
    let mut it = c.iter().rev();
    // accumulators for sum
    let mut y0 = if c.len() & 1 == 1 {
        *it.next().unwrap()
    } else {
        0.0
    };
    let mut y1 = 0.0;
    // Now n is even
    while let (Some(&c1), Some(&c0)) = (it.next(), it.next()) {
        // Unroll loop x 2, so accumulators return to their original role
        y1 = ar * y0 - y1 + c1;
        y0 = ar * y1 - y0 + c0;
    }
    (y0, y1)
}

/// Solve `k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0` for positive root k.
///
/// This solution is adapted from Geocentric::Reverse.
fn astroid(x: f64, y: f64) -> f64 {
    let p = sq(x);
    let q = sq(y);
    let r = (p + q - 1.0) / 6.0;
    if !(q == 0.0 && r <= 0.0) {
        // Avoid possible division by zero when r = 0 by multiplying equations
        // for s and t by r^3 and r, resp.
        let s = p * q / 4.0; // S = r^3 * s
        let r2 = sq(r);
        let r3 = r * r2;
        // The discriminant of the quadratic equation for T3.  This is zero on
        // the evolute curve p^(1/3)+q^(1/3) = 1
        let disc = s * (s + 2.0 * r3);
        let mut u = r;
        if disc >= 0.0 {
            let mut t3 = s + r3;
            // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
            // loss of precision due to cancellation.  The result is unchanged
            // because of the way the T is used in definition of u.
            t3 += if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() }; // T3 = (r * t)^3
            // N.B. cbrt always returns the f64 root.  cbrt(-8) = -2.
            let t = t3.cbrt(); // T = r * t
            // T can be zero; but then r2 / T -> 0.
            u += t + if t != 0.0 { r2 / t } else { 0.0 };
        } else {
            // T is complex, but the way u is defined the result is f64.
            let ang = (-disc).sqrt().atan2(-(s + r3));
            // There are three possible cube roots.  We choose the root which
            // avoids cancellation.  Note that disc < 0 implies that r < 0.
            u += 2.0 * r * (ang / 3.0).cos();
        }
        let v = (sq(u) + q).sqrt(); // guaranteed positive
        // Avoid loss of accuracy when u < 0.
        let uv = if u < 0.0 { q / (v - u) } else { u + v }; // u+v, guaranteed positive
        let w = (uv - q) / (2.0 * v); // positive?
        // Rearrange expression for k to avoid loss of accuracy due to
        // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
        uv / ((uv + sq(w)).sqrt() + w) // guaranteed positive
    } else {
        // q == 0 && r <= 0
        // y = 0 with |x| <= 1.  Handle this case directly.
        // for y small, positive root is k = abs(y)/sqrt(1-x^2)
        0.0
    }
}

/// Evaluate `coeff[o..=o+m] / coeff[o+m+1]` as a polynomial
#[inline]
fn polyval_ratio(coeff: &[f64], o: usize, m: usize, x: f64) -> f64 {
    polyval(&coeff[o..=o + m], x) / coeff[o + m + 1]
}

/// The scale factor A1-1 = mean value of (d/dsigma)I1 - 1
pub(crate) fn a1m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [
        // (1-eps)*A1-1, polynomial in eps2 of order 3
        1., 4., 64., 0., 256.,
    ];
    let m = NA1 / 2;
    let t = polyval_ratio(&COEFF, 0, m, sq(eps));
    (t + eps) / (1.0 - eps)
}

/// Evaluate the coefficients of a Fourier expansion as polynomials in eps^2
fn fourier_coeffs(coeff: &[f64], eps: f64) -> [f64; NC] {
    let mut c = [0.0; NC];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for (l, cl) in c.iter_mut().enumerate().skip(1) {
        let m = (NC - 1 - l) / 2; // order of polynomial in eps^2
        *cl = d * polyval(&coeff[o..=o + m], eps2) / coeff[o + m + 1];
        o += m + 2;
        d *= eps;
    }
    c
}

/// The coefficients C1[l] in the Fourier expansion of B1
pub(crate) fn c1f(eps: f64) -> [f64; NC] {
    const COEFF: [f64; 18] = [
        // C1[1]/eps^1, polynomial in eps2 of order 2
        -1., 6., -16., 32., //
        // C1[2]/eps^2, polynomial in eps2 of order 2
        -9., 64., -128., 2048., //
        // C1[3]/eps^3, polynomial in eps2 of order 1
        9., -16., 768., //
        // C1[4]/eps^4, polynomial in eps2 of order 1
        3., -5., 512., //
        // C1[5]/eps^5, polynomial in eps2 of order 0
        -7., 1280., //
        // C1[6]/eps^6, polynomial in eps2 of order 0
        -7., 2048.,
    ];
    fourier_coeffs(&COEFF, eps)
}

/// The coefficients C1p[l] in the Fourier expansion of B1p
pub(crate) fn c1pf(eps: f64) -> [f64; NC] {
    const COEFF: [f64; 18] = [
        // C1p[1]/eps^1, polynomial in eps2 of order 2
        205., -432., 768., 1536., //
        // C1p[2]/eps^2, polynomial in eps2 of order 2
        4005., -4736., 3840., 12288., //
        // C1p[3]/eps^3, polynomial in eps2 of order 1
        -225., 116., 384., //
        // C1p[4]/eps^4, polynomial in eps2 of order 1
        -7173., 2695., 7680., //
        // C1p[5]/eps^5, polynomial in eps2 of order 0
        3467., 7680., //
        // C1p[6]/eps^6, polynomial in eps2 of order 0
        38081., 61440.,
    ];
    fourier_coeffs(&COEFF, eps)
}

/// The scale factor A2-1 = mean value of (d/dsigma)I2 - 1
pub(crate) fn a2m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [
        // (eps+1)*A2-1, polynomial in eps2 of order 3
        -11., -28., -192., 0., 256.,
    ];
    let m = NA2 / 2;
    let t = polyval_ratio(&COEFF, 0, m, sq(eps));
    (t - eps) / (1.0 + eps)
}

/// The coefficients C2[l] in the Fourier expansion of B2
pub(crate) fn c2f(eps: f64) -> [f64; NC] {
    const COEFF: [f64; 18] = [
        // C2[1]/eps^1, polynomial in eps2 of order 2
        1., 2., 16., 32., //
        // C2[2]/eps^2, polynomial in eps2 of order 2
        35., 64., 384., 2048., //
        // C2[3]/eps^3, polynomial in eps2 of order 1
        15., 80., 768., //
        // C2[4]/eps^4, polynomial in eps2 of order 1
        7., 35., 512., //
        // C2[5]/eps^5, polynomial in eps2 of order 0
        63., 1280., //
        // C2[6]/eps^6, polynomial in eps2 of order 0
        77., 2048.,
    ];
    fourier_coeffs(&COEFF, eps)
}

/// The scale factor A3 = mean value of (d/dsigma)I3
fn a3_coeff(n: f64) -> [f64; NA3X] {
    const COEFF: [f64; 18] = [
        // A3, coeff of eps^5, polynomial in n of order 0
        -3., 128., //
        // A3, coeff of eps^4, polynomial in n of order 1
        -2., -3., 64., //
        // A3, coeff of eps^3, polynomial in n of order 2
        -1., -3., -1., 16., //
        // A3, coeff of eps^2, polynomial in n of order 2
        3., -1., -2., 8., //
        // A3, coeff of eps^1, polynomial in n of order 1
        1., -1., 2., //
        // A3, coeff of eps^0, polynomial in n of order 0
        1., 1.,
    ];
    let mut a3x = [0.0; NA3X];
    let mut o = 0;
    for (k, j) in (0..NA3).rev().enumerate() {
        // coeff of eps^j
        let m = (NA3 - j - 1).min(j); // order of polynomial in n
        a3x[k] = polyval_ratio(&COEFF, o, m, n);
        o += m + 2;
    }
    a3x
}

/// The coefficients C3[l] in the Fourier expansion of B3
fn c3_coeff(n: f64) -> [f64; NC3X] {
    const COEFF: [f64; 45] = [
        // C3[1], coeff of eps^5, polynomial in n of order 0
        3., 128., //
        // C3[1], coeff of eps^4, polynomial in n of order 1
        2., 5., 128., //
        // C3[1], coeff of eps^3, polynomial in n of order 2
        -1., 3., 3., 64., //
        // C3[1], coeff of eps^2, polynomial in n of order 2
        -1., 0., 1., 8., //
        // C3[1], coeff of eps^1, polynomial in n of order 1
        -1., 1., 4., //
        // C3[2], coeff of eps^5, polynomial in n of order 0
        5., 256., //
        // C3[2], coeff of eps^4, polynomial in n of order 1
        1., 3., 128., //
        // C3[2], coeff of eps^3, polynomial in n of order 2
        -3., -2., 3., 64., //
        // C3[2], coeff of eps^2, polynomial in n of order 2
        1., -3., 2., 32., //
        // C3[3], coeff of eps^5, polynomial in n of order 0
        7., 512., //
        // C3[3], coeff of eps^4, polynomial in n of order 1
        -10., 9., 384., //
        // C3[3], coeff of eps^3, polynomial in n of order 2
        5., -9., 5., 192., //
        // C3[4], coeff of eps^5, polynomial in n of order 0
        7., 512., //
        // C3[4], coeff of eps^4, polynomial in n of order 1
        -14., 7., 512., //
        // C3[5], coeff of eps^5, polynomial in n of order 0
        21., 2560.,
    ];
    let mut c3x = [0.0; NC3X];
    let (mut o, mut k) = (0, 0);
    for l in 1..NC3 {
        // l is index of C3[l]
        for j in (l..NC3).rev() {
            // coeff of eps^j
            let m = (NC3 - j - 1).min(j); // order of polynomial in n
            c3x[k] = polyval_ratio(&COEFF, o, m, n);
            k += 1;
            o += m + 2;
        }
    }
    c3x
}

/// The coefficients C4[l] in the Fourier expansion of I4
fn c4_coeff(n: f64) -> [f64; NC4X] {
    const COEFF: [f64; 77] = [
        // C4[0], coeff of eps^5, polynomial in n of order 0
        97., 15015., //
        // C4[0], coeff of eps^4, polynomial in n of order 1
        1088., 156., 45045., //
        // C4[0], coeff of eps^3, polynomial in n of order 2
        -224., -4784., 1573., 45045., //
        // C4[0], coeff of eps^2, polynomial in n of order 3
        -10656., 14144., -4576., -858., 45045., //
        // C4[0], coeff of eps^1, polynomial in n of order 4
        64., 624., -4576., 6864., -3003., 15015., //
        // C4[0], coeff of eps^0, polynomial in n of order 5
        100., 208., 572., 3432., -12012., 30030., 45045., //
        // C4[1], coeff of eps^5, polynomial in n of order 0
        1., 9009., //
        // C4[1], coeff of eps^4, polynomial in n of order 1
        -2944., 468., 135135., //
        // C4[1], coeff of eps^3, polynomial in n of order 2
        5792., 1040., -1287., 135135., //
        // C4[1], coeff of eps^2, polynomial in n of order 3
        5952., -11648., 9152., -2574., 135135., //
        // C4[1], coeff of eps^1, polynomial in n of order 4
        -64., -624., 4576., -6864., 3003., 135135., //
        // C4[2], coeff of eps^5, polynomial in n of order 0
        8., 10725., //
        // C4[2], coeff of eps^4, polynomial in n of order 1
        1856., -936., 225225., //
        // C4[2], coeff of eps^3, polynomial in n of order 2
        -8448., 4992., -1144., 225225., //
        // C4[2], coeff of eps^2, polynomial in n of order 3
        -1440., 4160., -4576., 1716., 225225., //
        // C4[3], coeff of eps^5, polynomial in n of order 0
        -136., 63063., //
        // C4[3], coeff of eps^4, polynomial in n of order 1
        1024., -208., 105105., //
        // C4[3], coeff of eps^3, polynomial in n of order 2
        3584., -3328., 1144., 315315., //
        // C4[4], coeff of eps^5, polynomial in n of order 0
        -128., 135135., //
        // C4[4], coeff of eps^4, polynomial in n of order 1
        -2560., 832., 405405., //
        // C4[5], coeff of eps^5, polynomial in n of order 0
        128., 99099.,
    ];
    let mut c4x = [0.0; NC4X];
    let (mut o, mut k) = (0, 0);
    for l in 0..NC4 {
        // l is index of C4[l]
        for j in (l..NC4).rev() {
            // coeff of eps^j
            let m = NC4 - j - 1; // order of polynomial in n
            c4x[k] = polyval_ratio(&COEFF, o, m, n);
            k += 1;
            o += m + 2;
        }
    }
    c4x
}
