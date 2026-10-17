//! Geodesic line
use crate::geodesic::*;
use crate::math::*;

/// A single geodesic, starting at point 1 with a given azimuth.
///
/// A `GeodesicLine` allows the efficient computation of many points
/// along the same geodesic.
///
/// ```rust
/// // Compute way points between JFK and Singapore Changi Airport,
/// // evenly spaced in arc length.
/// use proj4rs_geodesic::{Flags, Geodesic, Mask};
/// let g = Geodesic::wgs84();
/// let l = g.inverse_line(40.64, -73.78, 1.36, 103.99, Mask::LATITUDE | Mask::LONGITUDE);
/// for i in 0..=100 {
///     let p = l.gen_position(Flags::ARC_MODE, i as f64 * l.a13() * 0.01, Mask::LATITUDE | Mask::LONGITUDE);
///     println!("{:.5} {:.5}", p.lat2, p.lon2);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct GeodesicLine {
    lat1: f64,
    lon1: f64,
    azi1: f64,
    a: f64,
    f: f64,
    salp1: f64,
    calp1: f64,
    a13: f64,
    s13: f64,
    b: f64,
    c2: f64,
    f1: f64,
    salp0: f64,
    calp0: f64,
    k2: f64,
    ssig1: f64,
    csig1: f64,
    dn1: f64,
    stau1: f64,
    ctau1: f64,
    somg1: f64,
    comg1: f64,
    a1m1: f64,
    a2m1: f64,
    a3c: f64,
    b11: f64,
    b21: f64,
    b31: f64,
    a4: f64,
    b41: f64,
    c1a: [f64; NC],
    c1pa: [f64; NC],
    c2a: [f64; NC],
    c3a: [f64; NC3],
    c4a: [f64; NC4],
    caps: Mask,
}

impl GeodesicLine {
    /// Create a geodesic line starting at (`lat1`, `lon1`) with azimuth `azi1`
    /// (all in degrees).
    ///
    /// `caps` specifies the capabilities the line should possess, i.e.
    /// which quantities can be returned by [`GeodesicLine::gen_position`]:
    ///  - [`Mask::LATITUDE`] for the latitude `lat2`; this is added automatically,
    ///  - [`Mask::LONGITUDE`] for the longitude `lon2`,
    ///  - [`Mask::AZIMUTH`] for the azimuth `azi2`; this is added automatically,
    ///  - [`Mask::DISTANCE`] for the distance `s12`,
    ///  - [`Mask::REDUCED_LENGTH`] for the reduced length `m12`,
    ///  - [`Mask::GEODESIC_SCALE`] for the geodesic scales `M12` and `M21`,
    ///  - [`Mask::AREA`] for the area `S12`,
    ///  - [`Mask::DISTANCE_IN`] permits the length of the geodesic to be given
    ///    in terms of `s12`; without this capability the length can only be
    ///    specified in terms of arc length.
    ///
    /// A value of `caps` = [`Mask::NONE`] is treated as
    /// `LATITUDE | LONGITUDE | AZIMUTH | DISTANCE_IN` (to support the solution
    /// of the "standard" direct problem).
    ///
    /// Point 3 is undefined (`s13` and `a13` are NaN).
    pub fn new(g: &Geodesic, lat1: f64, lon1: f64, azi1: f64, caps: Mask) -> Self {
        let azi1 = ang_normalize(azi1);
        // Guard against underflow in salp0
        let (salp1, calp1) = sincosd(ang_round(azi1));
        Self::new_int(g, lat1, lon1, azi1, salp1, calp1, caps)
    }

    pub(crate) fn new_int(
        g: &Geodesic,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        salp1: f64,
        calp1: f64,
        caps: Mask,
    ) -> Self {
        // If caps is 0 assume the standard direct calculation
        let caps = if caps == Mask::NONE {
            Mask::DISTANCE_IN | Mask::LONGITUDE
        } else {
            caps
        } |
            // always allow latitude and azimuth
            Mask::LATITUDE | Mask::AZIMUTH;

        let lat1 = lat_fix(lat1);

        let (sbet1, cbet1) = sincosd(ang_round(lat1));
        // Ensure cbet1 = +epsilon at poles
        let (sbet1, cbet1) = norm2(sbet1 * g.f1, cbet1);
        let cbet1 = cbet1.max(TINY);
        let dn1 = (1.0 + g.ep2 * sq(sbet1)).sqrt();

        // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
        let salp0 = salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
        // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
        // is slightly better (consider the case salp1 = 0).
        let calp0 = calp1.hypot(salp1 * sbet1);
        // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
        // sig = 0 is nearest northward crossing of equator.
        // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
        // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
        // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
        // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
        // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
        // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
        // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
        let somg1 = salp0 * sbet1;
        let comg1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let (ssig1, csig1) = norm2(sbet1, comg1); // sig1 in (-pi, pi]
        // norm2(somg1, comg1); -- don't need to normalize!

        let k2 = sq(calp0) * g.ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);

        let mut l = Self {
            lat1,
            lon1,
            azi1,
            a: g.a,
            f: g.f,
            salp1,
            calp1,
            a13: f64::NAN,
            s13: f64::NAN,
            b: g.b,
            c2: g.c2,
            f1: g.f1,
            salp0,
            calp0,
            k2,
            ssig1,
            csig1,
            dn1,
            stau1: 0.0,
            ctau1: 0.0,
            somg1,
            comg1,
            a1m1: 0.0,
            a2m1: 0.0,
            a3c: 0.0,
            b11: 0.0,
            b21: 0.0,
            b31: 0.0,
            a4: 0.0,
            b41: 0.0,
            c1a: [0.0; NC],
            c1pa: [0.0; NC],
            c2a: [0.0; NC],
            c3a: [0.0; NC3],
            c4a: [0.0; NC4],
            caps,
        };

        if caps.any(Mask::CAP_C1) {
            l.a1m1 = a1m1f(eps);
            l.c1a = c1f(eps);
            l.b11 = sin_series(ssig1, csig1, &l.c1a[1..=NC1]);
            let (s, c) = l.b11.sin_cos();
            // tau1 = sig1 + B11
            l.stau1 = ssig1 * c + csig1 * s;
            l.ctau1 = csig1 * c - ssig1 * s;
            // Not necessary because C1pa reverts C1a
            //    B11 = -sin_series(stau1, ctau1, C1pa); */
        }

        if caps.any(Mask::CAP_C1P) {
            l.c1pa = c1pf(eps);
        }

        if caps.any(Mask::CAP_C2) {
            l.a2m1 = a2m1f(eps);
            l.c2a = c2f(eps);
            l.b21 = sin_series(ssig1, csig1, &l.c2a[1..=NC2]);
        }

        if caps.any(Mask::CAP_C3) {
            l.c3a = g.c3f(eps);
            l.a3c = -l.f * salp0 * g.a3f(eps);
            l.b31 = sin_series(ssig1, csig1, &l.c3a[1..NC3]);
        }

        if caps.any(Mask::CAP_C4) {
            l.c4a = g.c4f(eps);
            // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
            l.a4 = sq(l.a) * calp0 * salp0 * g.e2;
            l.b41 = cos_series(ssig1, csig1, &l.c4a);
        }

        l
    }

    /// Latitude of point 1 \[degrees\]
    pub fn lat1(&self) -> f64 {
        self.lat1
    }

    /// Longitude of point 1 \[degrees\]
    pub fn lon1(&self) -> f64 {
        self.lon1
    }

    /// Azimuth at point 1 \[degrees\]
    pub fn azi1(&self) -> f64 {
        self.azi1
    }

    /// Equatorial radius of the ellipsoid \[meters\]
    pub fn a(&self) -> f64 {
        self.a
    }

    /// Flattening of the ellipsoid
    pub fn f(&self) -> f64 {
        self.f
    }

    /// Arc length to the reference point 3 \[degrees\]
    pub fn a13(&self) -> f64 {
        self.a13
    }

    /// Distance to the reference point 3 \[meters\]
    pub fn s13(&self) -> f64 {
        self.s13
    }

    /// The capabilities of the line
    pub fn caps(&self) -> Mask {
        self.caps
    }

    /// Compute the position at distance `s12` \[meters\] along the line.
    ///
    /// Returns `(lat2, lon2, azi2)` in degrees; `lon2` is NaN if the
    /// line was not created with [`Mask::LONGITUDE`] and all values are
    /// NaN if it was not created with [`Mask::DISTANCE_IN`].
    pub fn position(&self, s12: f64) -> (f64, f64, f64) {
        let p = self.gen_position(
            Flags::NONE,
            s12,
            Mask::LATITUDE | Mask::LONGITUDE | Mask::AZIMUTH,
        );
        (p.lat2, p.lon2, p.azi2)
    }

    /// The general position function.
    ///
    /// If `flags` contains [`Flags::ARC_MODE`], `s12_a12` is the arc length
    /// from point 1 \[degrees\], otherwise it is the distance \[meters\]
    /// and the line must have been created with [`Mask::DISTANCE_IN`].
    /// With [`Flags::LONG_UNROLL`] the longitude is "unrolled" so that
    /// `lon2 - lon1` indicates how many times and in what sense the geodesic
    /// encircles the ellipsoid.
    ///
    /// Requesting a value which the line is not capable of computing is not
    /// an error; the corresponding value is NaN.
    pub fn gen_position(&self, flags: Flags, s12_a12: f64, outmask: Mask) -> Position {
        let mut r = Position::default();
        let outmask = outmask & self.caps & Mask(Mask::OUT_ALL);
        let arcmode = flags.contains(Flags::ARC_MODE);
        if !(arcmode || self.caps.any(Mask::DISTANCE_IN.out())) {
            // Impossible distance calculation requested
            return r;
        }

        let mut b12 = 0.0;
        let mut ab1 = 0.0;
        let (mut sig12, mut ssig12, mut csig12);
        if arcmode {
            // Interpret s12_a12 as spherical arc length
            sig12 = s12_a12 * DEGREE;
            (ssig12, csig12) = sincosd(s12_a12);
        } else {
            // Interpret s12_a12 as distance
            let tau12 = s12_a12 / (self.b * (1.0 + self.a1m1));
            let (s, c) = tau12.sin_cos();
            // tau2 = tau1 + tau12
            b12 = -sin_series(
                self.stau1 * c + self.ctau1 * s,
                self.ctau1 * c - self.stau1 * s,
                &self.c1pa[1..=NC1P],
            );
            sig12 = tau12 - (b12 - self.b11);
            (ssig12, csig12) = sig12.sin_cos();
            if self.f.abs() > 0.01 {
                // Reverted distance series is inaccurate for |f| > 1/100, so
                // correct sig12 with 1 Newton iteration.  The following table
                // shows the approximate maximum error for a = WGS_a() and
                // various f relative to GeodesicExact.
                //     erri = the error in the inverse solution (nm)
                //     errd = the error in the direct solution (series only) (nm)
                //     errda = the error in the direct solution
                //             (series + 1 Newton) (nm)
                //
                //       f     erri  errd errda
                //     -1/5    12e6 1.2e9  69e6
                //     -1/10  123e3  12e6 765e3
                //     -1/20   1110 108e3  7155
                //     -1/50  18.63 200.9 27.12
                //     -1/100 18.63 23.78 23.37
                //     -1/150 18.63 21.05 20.26
                //      1/150 22.35 24.73 25.83
                //      1/100 22.35 25.03 25.31
                //      1/50  29.80 231.9 30.44
                //      1/20   5376 146e3  10e3
                //      1/10  829e3  22e6 1.5e6
                //      1/5   157e6 3.8e9 280e6
                let ssig2 = self.ssig1 * csig12 + self.csig1 * ssig12;
                let csig2 = self.csig1 * csig12 - self.ssig1 * ssig12;
                b12 = sin_series(ssig2, csig2, &self.c1a[1..=NC1]);
                let serr = (1.0 + self.a1m1) * (sig12 + (b12 - self.b11)) - s12_a12 / self.b;
                sig12 -= serr / (1.0 + self.k2 * sq(ssig2)).sqrt();
                (ssig12, csig12) = sig12.sin_cos();
                // Update B12 below
            }
        }

        // sig2 = sig1 + sig12
        let ssig2 = self.ssig1 * csig12 + self.csig1 * ssig12;
        let mut csig2 = self.csig1 * csig12 - self.ssig1 * ssig12;
        let dn2 = (1.0 + self.k2 * sq(ssig2)).sqrt();
        if outmask.any((Mask::DISTANCE | Mask::REDUCED_LENGTH | Mask::GEODESIC_SCALE).out()) {
            if arcmode || self.f.abs() > 0.01 {
                b12 = sin_series(ssig2, csig2, &self.c1a[1..=NC1]);
            }
            ab1 = (1.0 + self.a1m1) * (b12 - self.b11);
        }
        // sin(bet2) = cos(alp0) * sin(sig2)
        let sbet2 = self.calp0 * ssig2;
        // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
        let mut cbet2 = self.salp0.hypot(self.calp0 * csig2);
        if cbet2 == 0.0 {
            // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
            cbet2 = TINY;
            csig2 = TINY;
        }
        // tan(alp0) = cos(sig2)*tan(alp2)
        let salp2 = self.salp0;
        let calp2 = self.calp0 * csig2; // No need to normalize

        if outmask.any(Mask::DISTANCE.out()) {
            r.s12 = if arcmode {
                self.b * ((1.0 + self.a1m1) * sig12 + ab1)
            } else {
                s12_a12
            };
        }

        if outmask.any(Mask::LONGITUDE.out()) {
            let e = 1.0_f64.copysign(self.salp0); // east or west going?
            // tan(omg2) = sin(alp0) * tan(sig2)
            let somg2 = self.salp0 * ssig2;
            let comg2 = csig2; // No need to normalize
            let unroll = flags.contains(Flags::LONG_UNROLL);
            // omg12 = omg2 - omg1
            let omg12 = if unroll {
                e * (sig12 - (ssig2.atan2(csig2) - self.ssig1.atan2(self.csig1))
                    + ((e * somg2).atan2(comg2) - (e * self.somg1).atan2(self.comg1)))
            } else {
                (somg2 * self.comg1 - comg2 * self.somg1)
                    .atan2(comg2 * self.comg1 + somg2 * self.somg1)
            };
            let lam12 = omg12
                + self.a3c * (sig12 + (sin_series(ssig2, csig2, &self.c3a[1..NC3]) - self.b31));
            let lon12 = lam12 / DEGREE;
            r.lon2 = if unroll {
                self.lon1 + lon12
            } else {
                ang_normalize(ang_normalize(self.lon1) + ang_normalize(lon12))
            };
        }

        if outmask.any(Mask::LATITUDE.out()) {
            r.lat2 = atan2d(sbet2, self.f1 * cbet2);
        }

        if outmask.any(Mask::AZIMUTH.out()) {
            r.azi2 = atan2d(salp2, calp2);
        }

        if outmask.any((Mask::REDUCED_LENGTH | Mask::GEODESIC_SCALE).out()) {
            let b22 = sin_series(ssig2, csig2, &self.c2a[1..=NC2]);
            let ab2 = (1.0 + self.a2m1) * (b22 - self.b21);
            let j12 = (self.a1m1 - self.a2m1) * sig12 + (ab1 - ab2);
            if outmask.any(Mask::REDUCED_LENGTH.out()) {
                // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
                // accurate cancellation in the case of coincident points.
                r.m12 = self.b
                    * ((dn2 * (self.csig1 * ssig2) - self.dn1 * (self.ssig1 * csig2))
                        - self.csig1 * csig2 * j12);
            }
            if outmask.any(Mask::GEODESIC_SCALE.out()) {
                let t = self.k2 * (ssig2 - self.ssig1) * (ssig2 + self.ssig1) / (self.dn1 + dn2);
                r.scale12 = csig12 + (t * ssig2 - csig2 * j12) * self.ssig1 / self.dn1;
                r.scale21 = csig12 - (t * self.ssig1 - self.csig1 * j12) * ssig2 / dn2;
            }
        }

        if outmask.any(Mask::AREA.out()) {
            let b42 = cos_series(ssig2, csig2, &self.c4a);
            let (salp12, calp12) = if self.calp0 == 0.0 || self.salp0 == 0.0 {
                // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                (
                    salp2 * self.calp1 - calp2 * self.salp1,
                    calp2 * self.calp1 + salp2 * self.salp1,
                )
            } else {
                // tan(alp) = tan(alp0) * sec(sig)
                // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
                // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
                // If csig12 > 0, write
                //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
                // else
                //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
                // No need to normalize
                (
                    self.calp0
                        * self.salp0
                        * (if csig12 <= 0.0 {
                            self.csig1 * (1.0 - csig12) + ssig12 * self.ssig1
                        } else {
                            ssig12 * (self.csig1 * ssig12 / (1.0 + csig12) + self.ssig1)
                        }),
                    sq(self.salp0) + sq(self.calp0) * self.csig1 * csig2,
                )
            };
            r.area12 = self.c2 * salp12.atan2(calp12) + self.a4 * (b42 - self.b41);
        }

        r.a12 = if arcmode { s12_a12 } else { sig12 / DEGREE };
        r
    }

    /// Specify position of point 3 in terms of distance \[meters\].
    ///
    /// This is only useful if the line has been created with
    /// [`Mask::DISTANCE_IN`].
    pub fn set_distance(&mut self, s13: f64) {
        self.s13 = s13;
        self.a13 = self.gen_position(Flags::NONE, s13, Mask::NONE).a12;
    }

    /// Specify position of point 3 in terms of arc length \[degrees\].
    ///
    /// The distance `s13` is only set if the line has been created
    /// with [`Mask::DISTANCE`].
    pub fn set_arc(&mut self, a13: f64) {
        self.a13 = a13;
        self.s13 = self.gen_position(Flags::ARC_MODE, a13, Mask::DISTANCE).s12;
    }

    /// Specify position of point 3 in terms of either distance \[meters\]
    /// or arc length \[degrees\] depending on [`Flags::ARC_MODE`].
    pub fn gen_set_distance(&mut self, flags: Flags, s13_a13: f64) {
        if flags.contains(Flags::ARC_MODE) {
            self.set_arc(s13_a13)
        } else {
            self.set_distance(s13_a13)
        }
    }
}
