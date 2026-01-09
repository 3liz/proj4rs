//! Rust interface for geodesic calculation from [GeographicLib](https://geographiclib.sourceforge.io/html/)
//! This is the original implementation in C with a Rust interface and inspired
//! from the [geographiclib rust project](https://github.com/savage13/geographiclib/tree/master)
//!
//! This is only a part of the geographiclib implementation needed by proj4rs
//! for implementing some projections.
//!
//! Example
//!
//! ```rust
//! use proj4rs_geodesic::Geodesic;
//! let g = Geodesic::wgs84();
//! let (lat1, lon1) = (37.87622, -122.23558); // Berkeley, California
//! let (lat2, lon2) = (-9.4047, 147.1597);    // Port Moresby, New Guinea
//! let (d_m, az1, az2) = g.inverse(lat1, lon1, lat2, lon2);
//!
//! assert_eq!(d_m, 10700471.955233702);  // Distance in meters
//! assert_eq!(az1, -96.91639942294974);  // Azimuth at (lat1, lon1)
//! assert_eq!(az2, -127.32548874543627); // Azimuth at (lat2, lon2)
//! ```
//!
//! Note: this is an alternative to the Vincenty distance calculation.
//!

use std::sync;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
static GEOD_INIT: std::cell::OnceCell<bool> = std::cell::OnceCell::new();

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
static GEOD_INIT: sync::OnceLock<bool> = sync::OnceLock::new();

/// Ellipsoid on which Geodesic Calculations are computed
#[repr(C)]
#[derive(Clone)]
pub struct Geodesic {
    /// Semi-major axis
    a: f64,
    /// Flattening
    f: f64,
    f1: f64,
    e2: f64,
    ep2: f64,
    n: f64,
    b: f64,
    c2: f64,
    etol2: f64,
    a3x: [f64; 6],
    c3x: [f64; 15],
    c4x: [f64; 21],
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

#[link(name = "geodesic", kind = "static")]
unsafe extern "C" {
    fn Init();
    fn geod_init(g: *mut Geodesic, a: f64, f: f64) -> bool;
    fn geod_inverse(
        g: *const Geodesic,
        lat1: f64,
        lon1: f64,
        lat2: f64,
        lon2: f64,
        ps12: *mut f64,
        pazi1: *mut f64,
        pazi2: *mut f64,
    );
    fn geod_direct(
        g: *const Geodesic,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        s12: f64,
        plat2: *mut f64,
        plon2: *mut f64,
        pazi2: *mut f64,
    );
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
        GEOD_INIT.get_or_init(|| {
            unsafe {
                Init();
            }
            true
        });
        unsafe {
            let mut g = std::mem::MaybeUninit::<Geodesic>::uninit();
            if !geod_init(g.as_mut_ptr(), a, f) {
                panic!("geodesic is not initialized");
            }
            g.assume_init()
        }
    }

    #[allow(non_upper_case_globals)]
    pub fn wgs84() -> Self {
        const a: f64 = 6_378_137.0;
        const f: f64 = 1.0 / 298.257_223_563; /* WGS84 */
        Self::new(a, f)
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
        let mut s12 = 0.0;
        let mut azi1 = 0.0;
        let mut azi2 = 0.0;
        unsafe {
            geod_inverse(
                self as *const Geodesic,
                lat1,
                lon1,
                lat2,
                lon2,
                &mut s12 as *mut f64,
                &mut azi1 as *mut f64,
                &mut azi2 as *mut f64,
            )
        };
        (s12, azi1, azi2)
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
        let mut lat2 = 0.0;
        let mut lon2 = 0.0;
        let mut azi2 = 0.0;
        unsafe {
            geod_direct(
                self as *const Geodesic,
                lat1,
                lon1,
                azi1,
                s12,
                &mut lat2 as &mut f64,
                &mut lon2 as &mut f64,
                &mut azi2 as &mut f64,
            )
        };
        (lat2, lon2, azi2)
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn dist_az_test() {
        struct TestCase {
            pub lat1: f64,
            pub lon1: f64,
            pub azi1: f64,
            pub lat2: f64,
            pub lon2: f64,
            pub azi2: f64,
            pub s12: f64,
            //pub a12: f64,
            //pub m12: f64,
            //pub mm12: f64, // M12
            //pub mm21: f64, // M21
            //pub ss12: f64, // S12
        }
        impl TestCase {
            fn vec(v: &[f64]) -> Self {
                Self {
                    lat1: v[0],
                    lon1: v[1],
                    azi1: v[2],
                    lat2: v[3],
                    lon2: v[4],
                    azi2: v[5],
                    s12: v[6],
                    //a12: v[7],
                    //m12: v[8],
                    //mm12: v[9],
                    //mm21: v[10],
                    //ss12: v[11],
                }
            }
        }

        let testcases = [
            TestCase::vec(&[
                35.60777,
                -139.44815,
                111.098748429560326,
                -11.17491,
                -69.95921,
                129.289270889708762,
                8935244.5604818305,
                80.50729714281974,
                6273170.2055303837,
                0.16606318447386067,
                0.16479116945612937,
                12841384694976.432,
            ]),
            TestCase::vec(&[
                55.52454,
                106.05087,
                22.020059880982801,
                77.03196,
                197.18234,
                109.112041110671519,
                4105086.1713924406,
                36.892740690445894,
                3828869.3344387607,
                0.80076349608092607,
                0.80101006984201008,
                61674961290615.615,
            ]),
            TestCase::vec(&[
                -21.97856,
                142.59065,
                -32.44456876433189,
                41.84138,
                98.56635,
                -41.84359951440466,
                8394328.894657671,
                75.62930491011522,
                6161154.5773110616,
                0.24816339233950381,
                0.24930251203627892,
                -6637997720646.717,
            ]),
            TestCase::vec(&[
                -66.99028,
                112.2363,
                173.73491240878403,
                -12.70631,
                285.90344,
                2.512956620913668,
                11150344.2312080241,
                100.278634181155759,
                6289939.5670446687,
                -0.17199490274700385,
                -0.17722569526345708,
                -121287239862139.744,
            ]),
            TestCase::vec(&[
                -17.42761,
                173.34268,
                -159.033557661192928,
                -15.84784,
                5.93557,
                -20.787484651536988,
                16076603.1631180673,
                144.640108810286253,
                3732902.1583877189,
                -0.81273638700070476,
                -0.81299800519154474,
                97825992354058.708,
            ]),
            TestCase::vec(&[
                32.84994,
                48.28919,
                150.492927788121982,
                -56.28556,
                202.29132,
                48.113449399816759,
                16727068.9438164461,
                150.565799985466607,
                3147838.1910180939,
                -0.87334918086923126,
                -0.86505036767110637,
                -72445258525585.010,
            ]),
            TestCase::vec(&[
                6.96833,
                52.74123,
                92.581585386317712,
                -7.39675,
                206.17291,
                90.721692165923907,
                17102477.2496958388,
                154.147366239113561,
                2772035.6169917581,
                -0.89991282520302447,
                -0.89986892177110739,
                -1311796973197.995,
            ]),
            TestCase::vec(&[
                -50.56724,
                -16.30485,
                -105.439679907590164,
                -33.56571,
                -94.97412,
                -47.348547835650331,
                6455670.5118668696,
                58.083719495371259,
                5409150.7979815838,
                0.53053508035997263,
                0.52988722644436602,
                41071447902810.047,
            ]),
            TestCase::vec(&[
                -58.93002,
                -8.90775,
                140.965397902500679,
                -8.91104,
                133.13503,
                19.255429433416599,
                11756066.0219864627,
                105.755691241406877,
                6151101.2270708536,
                -0.26548622269867183,
                -0.27068483874510741,
                -86143460552774.735,
            ]),
            TestCase::vec(&[
                -68.82867,
                -74.28391,
                93.774347763114881,
                -50.63005,
                -8.36685,
                34.65564085411343,
                3956936.926063544,
                35.572254987389284,
                3708890.9544062657,
                0.81443963736383502,
                0.81420859815358342,
                -41845309450093.787,
            ]),
            TestCase::vec(&[
                -10.62672,
                -32.0898,
                -86.426713286747751,
                5.883,
                -134.31681,
                -80.473780971034875,
                11470869.3864563009,
                103.387395634504061,
                6184411.6622659713,
                -0.23138683500430237,
                -0.23155097622286792,
                4198803992123.548,
            ]),
            TestCase::vec(&[
                -21.76221,
                166.90563,
                29.319421206936428,
                48.72884,
                213.97627,
                43.508671946410168,
                9098627.3986554915,
                81.963476716121964,
                6299240.9166992283,
                0.13965943368590333,
                0.14152969707656796,
                10024709850277.476,
            ]),
            TestCase::vec(&[
                -19.79938,
                -174.47484,
                71.167275780171533,
                -11.99349,
                -154.35109,
                65.589099775199228,
                2319004.8601169389,
                20.896611684802389,
                2267960.8703918325,
                0.93427001867125849,
                0.93424887135032789,
                -3935477535005.785,
            ]),
            TestCase::vec(&[
                -11.95887,
                -116.94513,
                92.712619830452549,
                4.57352,
                7.16501,
                78.64960934409585,
                13834722.5801401374,
                124.688684161089762,
                5228093.177931598,
                -0.56879356755666463,
                -0.56918731952397221,
                -9919582785894.853,
            ]),
            TestCase::vec(&[
                -87.85331,
                85.66836,
                -65.120313040242748,
                66.48646,
                16.09921,
                -4.888658719272296,
                17286615.3147144645,
                155.58592449699137,
                2635887.4729110181,
                -0.90697975771398578,
                -0.91095608883042767,
                42667211366919.534,
            ]),
            TestCase::vec(&[
                1.74708,
                128.32011,
                -101.584843631173858,
                -11.16617,
                11.87109,
                -86.325793296437476,
                12942901.1241347408,
                116.650512484301857,
                5682744.8413270572,
                -0.44857868222697644,
                -0.44824490340007729,
                10763055294345.653,
            ]),
            TestCase::vec(&[
                -25.72959,
                -144.90758,
                -153.647468693117198,
                -57.70581,
                -269.17879,
                -48.343983158876487,
                9413446.7452453107,
                84.664533838404295,
                6356176.6898881281,
                0.09492245755254703,
                0.09737058264766572,
                74515122850712.444,
            ]),
            TestCase::vec(&[
                -41.22777,
                122.32875,
                14.285113402275739,
                -7.57291,
                130.37946,
                10.805303085187369,
                3812686.035106021,
                34.34330804743883,
                3588703.8812128856,
                0.82605222593217889,
                0.82572158200920196,
                -2456961531057.857,
            ]),
            TestCase::vec(&[
                11.01307,
                138.25278,
                79.43682622782374,
                6.62726,
                247.05981,
                103.708090215522657,
                11911190.819018408,
                107.341669954114577,
                6070904.722786735,
                -0.29767608923657404,
                -0.29785143390252321,
                17121631423099.696,
            ]),
            TestCase::vec(&[
                -29.47124,
                95.14681,
                -163.779130441688382,
                -27.46601,
                -69.15955,
                -15.909335945554969,
                13487015.8381145492,
                121.294026715742277,
                5481428.9945736388,
                -0.51527225545373252,
                -0.51556587964721788,
                104679964020340.318,
            ]),
        ];
        let g = crate::Geodesic::wgs84();
        for t in &testcases {
            let (s, azi1, azi2) = g.inverse(t.lat1, t.lon1, t.lat2, t.lon2);
            assert!((s - t.s12).abs() < 1e-8, "{} {}", s, t.s12);
            assert!((azi1 - t.azi1).abs() < 1e-13, "{} {}", azi1, t.azi1);
            assert!((azi2 - t.azi2).abs() < 1e-13, "{} {}", azi2, t.azi2);
        }
        let (s, az1, az2) = g.inverse(0.0, 0.0, 0.0, 10.0);
        let s0 = 1113194.9079327357;
        assert!((s - s0).abs() < 1e-5, "{} {}", s, s0);
        assert_eq!(az1, 90.0);
        assert_eq!(az2, 90.0);
    }

    #[test]
    fn test_debug() {
        use crate::Geodesic;
        let g = Geodesic::new(6_378_145.0, 1.0 / 298.25);
        assert_eq!(
            format!("{}", g),
            "Geodesic { a: 6378145, f: 0.003352891869237217 }"
        );
    }
}
