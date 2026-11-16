//! Precision of the spherical projections close to the poles, the equator and the projection
//! centres. Expected values are from a 50-digit computation on the same input, on a sphere of
//! radius 6371000 m. The inverse is checked as a great-circle distance in metres, because close
//! to the poles the longitude is ill-conditioned.
//!
//! Each case also carries the output of PROJ 9.8.1 (pyproj) on the same input: forward (x, y) and
//! inverse (lon, lat). proj4rs must match it within 1e-6 m, except where PROJ itself is 1e-6 m or
//! more off the exact value (most of these points were chosen where PROJ loses precision); there,
//! proj4rs must be closer to the exact value than PROJ. See `check_proj()`.

use proj4rs::Proj;
use proj4rs::transform::transform;

const R: f64 = 6371000.;

/// Tolerance in metres against the PROJ 9.8.1 outputs
const TOL_PROJ: f64 = 1e-6;

/// Compares with the PROJ 9.8.1 output: `diff` is the distance in metres between proj4rs and PROJ,
/// `ours` and `theirs` the distances of proj4rs and PROJ to the exact value. Within TOL_PROJ of
/// PROJ passes; beyond it, PROJ must itself be at least TOL_PROJ off the exact value and proj4rs
/// closer to it than PROJ (a known PROJ precision loss).
fn check_proj(failures: &mut Vec<String>, what: String, diff: f64, ours: f64, theirs: f64) {
    if diff < TOL_PROJ || (theirs >= TOL_PROJ && ours < theirs) {
        return;
    }
    failures.push(format!(
        "{what}: {diff:e} m from PROJ 9.8.1 (proj4rs {ours:e} m, PROJ {theirs:e} m from the exact value)"
    ));
}

/// (def, lon, lat, exact x, exact y, PROJ x, PROJ y, PROJ inverse lon, PROJ inverse lat)
type Case = (&'static str, f64, f64, f64, f64, f64, f64, f64, f64);

#[rustfmt::skip]
const CASES: &[Case] = &[
    // stere sphere, polar: the inverse took the latitude from asin() near 1 (main: 90, 1.1 cm off)
    ("+proj=stere +lat_0=90 +R=6371000 +k_0=0.994 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.9999999, 0.005526387526143001, -0.009571983977594557, 0.005526387735470094, -0.009571984340159718, 29.999999999999996, 90.0),
    // ortho sphere, oblique: the inverse took the latitude from asin() near 1 (main: 1.1 cm off)
    ("+proj=ortho +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.9999999, 0.0019308795231307598, 4504977.295196278, 0.0019308796640102259, 4504977.295196277, 10.000000110531504, 90.0),
    // aeqd sphere, oblique: the inverse took the latitude from asin() near 1 (main: 1.1 cm off)
    ("+proj=aeqd +lat_0=-30 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -120.0, -89.9999999, -0.011644304801302343, -6671695.6042332705, -0.011644305650885879, -6671695.604233271, -120.0000016524417, -90.0),
    // laea sphere, oblique: the inverse took the latitude from asin() near 1 (main: 1.1 cm off)
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.9999999, 0.01042317987605592, 4876152.290059438, 0.010423180636544597, 4876152.290059438, 59.999998954543194, 90.0),
    // laea sphere, equatorial: within 1e-10 * R of the centre the inverse returned latitude 0
    // (main: 0.22 mm off)
    ("+proj=laea +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 1e-09, -2e-09, 0.00011119492664455874, -0.00022238985328911748, 0.00011119492664455875, -0.0002223898532891175, 9.999999999999999e-10, 0.0),
    // gnom sphere, polar: regression test for the polar branch (not changed; exact on main too)
    ("+proj=gnom +lat_0=-90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 45.0, -89.9999999, 0.007862668199598521, 0.007862668199598521, 0.007862668773268957, 0.007862668773268959, 45.0, -89.9999999),
    // gnom: the centre itself (regression test; main returned lat_0 for rh <= EPS_10 too)
    ("+proj=gnom +lat_0=90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 90.0, 0.0, 0.0, 0.0, -3.901112378683894e-10, 0.0, 90.0),
    // tpers: the inverse took the latitude from asin() near 1 (main: 1.1 cm off)
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.9999999, 0.0067512397377098715, 2777408.316916929, 0.0067512402302890495, 2777408.3169169286, 99.99999473497667, 90.0),
    // tpers: the centre itself (regression test; main returned lat_0 too)
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 45.0, 0.0, 0.0, 0.0, 0.0, 0.0, 45.0),
    // moll: near the poles the forward iteration stopped early (main: 591 m off) and the inverse
    // clamped sin(theta) (main: 4.8 mm off)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -179.0, -89.99999, -591.078391808469, -9009954.60097763, -591.223420892261, -9009954.600975225, -178.99998593096112, -89.99999004337388),
    // moll: near the poles the forward iteration stopped early (main: 33 km off)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.99, 33021.10324921167, 9009905.592261195, 33021.10323607112, 9009905.592261236, 99.9999999984145, 89.99000000001735),
    // tmerc sphere: close to the equator acos() lost half the digits (main: 1.1 cm off, both ways)
    ("+proj=tmerc +R=6371000 +k=0.999 +x_0=0 +y_0=0 +lon_0=0 +lat_0=0", 2.0, 1e-07, 222212.5945705513, 0.011115144217327467, 222212.59457055124, 0.0, 2.0000000000000058, 0.0),
    // laea sphere, oblique: the inverse gave longitude 0 instead of 90 (main: 9815 km off)
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 90.0, -10.0, 9473745.882611485, -1181205.6198054145, 9473745.882611483, -1181205.619805415, 0.0, -10.000000000000002),
    // eqearth: MAX_Y was too small, so the inverse at the pole was 0.74 m off on main
    ("+proj=eqearth +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 30.0, 90.0, 1702743.8797029483, 8392918.138591882, 1702743.8797029485, 8392918.13859188, 29.999999999999417, 89.99999338669778),
];

fn distance(l1: f64, p1: f64, l2: f64, p2: f64) -> f64 {
    let a = ((p2 - p1) / 2.).sin();
    let b = ((l2 - l1) / 2.).sin();
    2. * R * (a * a + p1.cos() * p2.cos() * b * b).sqrt().asin()
}

#[test]
fn precision_near_poles_and_centres() {
    let ll = Proj::from_proj_string("+proj=latlong +R=6371000").unwrap();
    let mut failures = Vec::new();
    for &(def, lon, lat, x, y, px, py, plon, plat) in CASES {
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (lon.to_radians(), lat.to_radians(), 0.);
        transform(&ll, &pj, &mut p).unwrap();
        let fwd = (p.0 - x).hypot(p.1 - y);
        if fwd.is_nan() || fwd >= 1e-6 {
            failures.push(format!(
                "{def} forward ({lon}, {lat}): ({}, {}) is {fwd:e} m off",
                p.0, p.1
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} forward ({lon}, {lat})"),
            (p.0 - px).hypot(p.1 - py),
            fwd,
            (px - x).hypot(py - y),
        );
        let mut q = (x, y, 0.);
        transform(&pj, &ll, &mut q).unwrap();
        let (lon, lat, plon, plat) = (
            lon.to_radians(),
            lat.to_radians(),
            plon.to_radians(),
            plat.to_radians(),
        );
        let inv = distance(q.0, q.1, lon, lat);
        if inv.is_nan() || inv >= 1e-6 {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {inv:e} m off",
                q.0.to_degrees(),
                q.1.to_degrees()
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} inverse ({x}, {y})"),
            distance(q.0, q.1, plon, plat),
            inv,
            distance(plon, plat, lon, lat),
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// gnom in equatorial and oblique mode, close to 90 degrees from the centre, where the projected
/// coordinates are very large and cos(z) computed as sqrt(1 - sin(z)^2) lost its digits: (x, y) is
/// the exact forward image rounded to double, and the expected (lon, lat) is the exact inverse of
/// those doubles, from a 50-digit computation.
/// (def, x, y, exact lon, exact lat, PROJ 9.8.1 inverse lon, PROJ 9.8.1 inverse lat)
#[rustfmt::skip]
const GNOM_FAR: &[(&str, f64, f64, f64, f64, f64, f64)] = &[
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -0.011119492664455874, -3650314329490979.0, -1e-07, -89.9999999, -90.0, -90.0),  // from (-1e-07, -89.9999999)
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 1123379.1940936304, 3706626311913009.5, 10.000000000000002, 89.9999999, 90.0, 90.0),  // from (10, 89.9999999)
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 36503141116197.78, 0.0, 89.99999, 0.0, 89.99998997043059),  // from (0, 89.99999)
    ("+proj=gnom +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 516232372521323.0, 6371000.0, 90.0, 1e-06, 89.99999989629086, 1.1037091348628663e-06),  // from (90, 1e-06)
    ("+proj=gnom +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, -3650314070120549.0, 0.0, -44.9999999, 0.0, -45.00000000000001),  // from (0, -44.9999999)
    ("+proj=gnom +lat_0=-30 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 3650314070120549.0, 0.0, 59.9999999, 0.0, 60.00000000000001),  // from (0, 59.9999999)
];

#[test]
fn gnom_inverse_far_from_the_centre() {
    let ll = Proj::from_proj_string("+proj=latlong +R=6371000").unwrap();
    let mut failures = Vec::new();
    for &(def, x, y, lon, lat, plon, plat) in GNOM_FAR {
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (x, y, 0.);
        transform(&pj, &ll, &mut p).unwrap();
        let (lon, lat, plon, plat) = (
            lon.to_radians(),
            lat.to_radians(),
            plon.to_radians(),
            plat.to_radians(),
        );
        let d = distance(p.0, p.1, lon, lat);
        if d.is_nan() || d >= 1e-6 {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {d:e} m off",
                p.0.to_degrees(),
                p.1.to_degrees()
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} inverse ({x}, {y})"),
            distance(p.0, p.1, plon, plat),
            d,
            distance(plon, plat, lon, lat),
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// moll: points on the map outline (longitude +/-180) must survive a round trip. The inverse
/// rejected longitudes that rounded to slightly more than pi there. How x rounds on the outline
/// depends on the platform's sin/cos, so x is also moved 1 and 2 ulps outwards.
#[test]
fn moll_outline_round_trip() {
    let ll = Proj::from_proj_string("+proj=latlong +R=6371000").unwrap();
    let pj = Proj::from_proj_string("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0").unwrap();
    // next double away from zero (the bits of an f64 grow with its magnitude, for either sign)
    let outwards = |v: f64| f64::from_bits(v.to_bits() + 1);
    let mut failures = Vec::new();
    for lon in [-180_f64, 180.] {
        for lat in [
            -89.9999999_f64,
            -75.,
            -45.,
            0.,
            30.,
            45.,
            61.4,
            75.,
            89.99,
            89.9999999,
        ] {
            let mut p = (lon.to_radians(), lat.to_radians(), 0.);
            transform(&ll, &pj, &mut p).unwrap();
            let y = p.1;
            for x in [p.0, outwards(p.0), outwards(outwards(p.0))] {
                let mut q = (x, y, 0.);
                if let Err(e) = transform(&pj, &ll, &mut q) {
                    failures.push(format!(
                        "({lon}, {lat}): inverse of ({x}, {y}) failed: {e:?}"
                    ));
                    continue;
                }
                let d = distance(q.0, q.1, lon.to_radians(), lat.to_radians());
                if d.is_nan() || d >= 1e-6 {
                    failures.push(format!(
                        "({lon}, {lat}): inverse of ({x}, {y}) gave ({}, {}), {d:e} m off",
                        q.0.to_degrees(),
                        q.1.to_degrees()
                    ));
                }
            }
            // 0.1 % outside the outline is still outside the projection
            #[cfg(not(feature = "proj4js-compat"))]
            {
                let xo = p.0 * 1.001;
                let mut q = (xo, y, 0.);
                if transform(&pj, &ll, &mut q).is_ok() {
                    failures.push(format!(
                        "({lon}, {lat}): x = {xo} outside the outline was accepted"
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

//
// Ellipsoidal projections (GRS80)
//

const A: f64 = 6378137.;

/// Distance in metres between two points given in radians, on a sphere of radius A: close enough
/// to the ellipsoidal distance for the sub-millimetre distances checked here.
fn distance_a(l1: f64, p1: f64, l2: f64, p2: f64) -> f64 {
    distance(l1, p1, l2, p2) * A / R
}

/// laea ellipsoidal, polar aspect: the forward took sqrt(qp - q), which cancels close to the pole,
/// and returned (0, 0) when qp - q < 1e-15, i.e. up to about 0.2 m from the pole. Expected values
/// are a * sqrt(qp -/+ q(phi)) * (sin(lam), -/+cos(lam)) from a 60-digit computation on the same
/// input; the errors of main (same in PROJ 9.8.1) are given for each case.
/// See `Case` for the fields.
#[rustfmt::skip]
const LAEA_POLAR: &[Case] = &[
    // main: (0, 0), 11 mm off
    ("+proj=laea +lat_0=90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.9999999, 0.005584699053961908, -0.009672982506443868, 0.0, 0.0, 29.999999999999996, 90.0),
    // main: (0, 0), 11 cm off
    ("+proj=laea +lat_0=90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.999999, 0.05584698948673728, -0.09672982324079392, 0.0, 0.0, 29.999999999999996, 89.9999991500401),
    // main: 0.45 mm off
    ("+proj=laea +lat_0=90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", -120.0, 89.99999, -0.9672982379680011, 0.5584698980774758, -0.9669057237710204, 0.5582432799001883, -119.99999999999999, 89.99999001523184),
    // main: (0, 0), 11 mm off
    ("+proj=laea +lat_0=-90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 45.0, -89.9999999, 0.007897957143885123, 0.007897957143885125, 0.0, 0.0, 45.0, -90.0),
    // main: (0, 0), 11 cm off
    ("+proj=laea +lat_0=-90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 150.0, -89.999999, 0.05584698948673728, -0.09672982324079392, 0.0, 0.0, 150.0, -89.9999991500401),
    // near the centre, main: (0, 0), 0.22 mm off
    ("+proj=laea +lat_0=-90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 1e-9, -89.999999998, 3.898855207191207e-15, 0.00022338794830466024, 0.0, 0.0, 1e-09, -90.0),
];

#[test]
fn laea_ellipsoidal_polar() {
    let ll = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
    let mut failures = Vec::new();
    for &(def, lon, lat, x, y, px, py, plon, plat) in LAEA_POLAR {
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (lon.to_radians(), lat.to_radians(), 0.);
        transform(&ll, &pj, &mut p).unwrap();
        let fwd = (p.0 - x).hypot(p.1 - y);
        if fwd.is_nan() || fwd >= 1e-9 {
            failures.push(format!(
                "{def} forward ({lon}, {lat}): ({}, {}) is {fwd:e} m off",
                p.0, p.1
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} forward ({lon}, {lat})"),
            (p.0 - px).hypot(p.1 - py),
            fwd,
            (px - x).hypot(py - y),
        );
        let mut q = (x, y, 0.);
        transform(&pj, &ll, &mut q).unwrap();
        let (lon, lat, plon, plat) = (
            lon.to_radians(),
            lat.to_radians(),
            plon.to_radians(),
            plat.to_radians(),
        );
        // 1e-8 m: one unit in the last place of the latitude is 1.4e-9 m close to the poles
        let inv = distance_a(q.0, q.1, lon, lat);
        if inv.is_nan() || inv >= 1e-8 {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {inv:e} m off",
                q.0.to_degrees(),
                q.1.to_degrees()
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} inverse ({x}, {y})"),
            distance_a(q.0, q.1, plon, plat),
            inv,
            distance_a(plon, plat, lon, lat),
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Ellipsoidal inverses close to the poles and the centres, checked as a round trip (the forward
/// is accurate at these points). laea took the authalic latitude from asin() of a value close to
/// +/-1 and returned the centre within 1e-10 * a of it; stere took the conformal latitude from
/// asin() close to +/-1 and kept that error through its iteration. The round-trip errors of main
/// (same in PROJ 9.8.1) are given for each case.
///
/// PROJ 9.8.1 outputs: its forward (x, y) of (lon, lat), which proj4rs must match within 1e-6 m
/// (the forward is accurate at these points), and its inverse of that (x, y), compared with the
/// inverse of proj4rs on the same (x, y), taking (lon, lat) as the exact value.
/// (def, lon, lat, PROJ x, PROJ y, PROJ inverse lon, PROJ inverse lat)
#[rustfmt::skip]
const ELLIPSOIDAL_ROUND_TRIP: &[(&str, f64, f64, f64, f64, f64, f64)] = &[
    // laea oblique, near the pole (main: 1.1 cm and 1.7 cm)
    ("+proj=laea +lat_0=45 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.9999999, 0.010480700226785903, 4886594.215564228, 60.00000319815302, 90.0),
    ("+proj=laea +lat_0=45 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.999999, 0.10480699614670957, 4886594.169175506, 59.999999709721834, 89.9999991500401),
    // laea equatorial, near the pole (main: 1.1 cm and 2.3 mm)
    ("+proj=laea +lat_0=0 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.9999999, 0.002746001465302188, 8999892.937928976, 9.999997616327905, 90.0),
    ("+proj=laea +lat_0=0 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.99999, 0.274600108350179, 8999892.168770744, 9.999999980847052, 89.99998997912043),
    // laea oblique, near the south pole (main: 8.3 cm)
    ("+proj=laea +lat_0=-30 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", -100.0, -89.9999999, -0.012719163656680315, -6376335.552622448, -100.00000016898294, -89.9999991500401),
    // laea near the centre: the inverse returned the centre (main: 0.25 mm and 0.24 mm)
    ("+proj=laea +lat_0=0 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 1e-9, -2e-9, 0.00011131949079327359, -0.00022114855163589525, 0.0, 0.0),
    ("+proj=laea +lat_0=45 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 1e-9, 44.999999998, 7.884683509736598e-05, -0.00022226363477656236, 0.0, 45.0),
    // stere oblique, near the pole (main: 1.1 cm and 2.2 cm)
    ("+proj=stere +lat_0=45 +k_0=1 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.9999999, 0.011367511663901266, 5300068.109808842, 59.999996238165004, 90.0),
    ("+proj=stere +lat_0=45 +k_0=1 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.999999, 0.1136751241918649, 5300068.050741513, 59.99999995552266, 89.99999880065565),
    // stere equatorial, near the pole (main: 1.1 cm and 1.7 cm)
    ("+proj=stere +lat_0=0 +k_0=1 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.9999999, 0.0038921112367458328, 12756273.97792674, 10.00000005177111, 90.0),
    ("+proj=stere +lat_0=0 +k_0=1 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.999999, 0.03892111447541261, 12756273.77926739, 10.000000041230782, 89.99999915193547),
    // stere oblique, near the south pole (main: 8.3 cm)
    ("+proj=stere +lat_0=-30 +k_0=0.9996 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", -100.0, -89.9999999, -0.014722070577685759, -7380439.899899116, -100.00000264890562, -89.99999915193547),
];

#[test]
fn ellipsoidal_round_trip_near_poles_and_centres() {
    let ll = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
    let mut failures = Vec::new();
    for &(def, lon, lat, px, py, plon, plat) in ELLIPSOIDAL_ROUND_TRIP {
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (lon.to_radians(), lat.to_radians(), 0.);
        transform(&ll, &pj, &mut p).unwrap();
        let (x, y) = (p.0, p.1);
        let fwd = (x - px).hypot(y - py);
        if fwd.is_nan() || fwd >= TOL_PROJ {
            failures.push(format!(
                "{def} forward ({lon}, {lat}): ({x}, {y}) is {fwd:e} m from PROJ 9.8.1 ({px}, {py})"
            ));
        }
        transform(&pj, &ll, &mut p).unwrap();
        let (lon, lat, plon, plat) = (
            lon.to_radians(),
            lat.to_radians(),
            plon.to_radians(),
            plat.to_radians(),
        );
        let d = distance_a(p.0, p.1, lon, lat);
        if d.is_nan() || d >= 1e-7 {
            failures.push(format!(
                "{def} ({}, {}) -> ({x}, {y}) -> ({}, {}) is {d:e} m off",
                lon.to_degrees(),
                lat.to_degrees(),
                p.0.to_degrees(),
                p.1.to_degrees()
            ));
        }
        // inverse of the same (x, y) as PROJ
        let mut q = (px, py, 0.);
        transform(&pj, &ll, &mut q).unwrap();
        check_proj(
            &mut failures,
            format!("{def} inverse ({px}, {py})"),
            distance_a(q.0, q.1, plon, plat),
            distance_a(q.0, q.1, lon, lat),
            distance_a(plon, plat, lon, lat),
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// aeqd close to the centre: the inverses (polar sphere, polar ellipsoid, geodesic) and the
/// geodesic forward returned the centre within 1e-10 (times the radius) of it. Exact values: R * c
/// on the sphere, a / sqrt(1 - e^2) * c (radius of curvature at the pole) in the ellipsoidal polar
/// aspect, and (N cos(phi0) dlam, M dphi) through the geodesic, where the neglected terms are below
/// 1e-15 m at these distances. The errors of main (same in PROJ 9.8.1) are given for each case.
/// See `Case` for the fields.
#[rustfmt::skip]
const AEQD_NEAR_CENTRE: &[Case] = &[
    // polar sphere, inverse (main: the pole, 0.11 mm and 0.22 mm off)
    ("+proj=aeqd +lat_0=90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.999999999, 5.559791181939671e-05, -9.629840806592932e-05, 5.559771676377777e-05, -9.629807021968702e-05, 0.0, 90.0),
    ("+proj=aeqd +lat_0=-90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -150.0, -89.999999998, -0.0001111949212600855, -0.00019259525316608884, -0.00011119472620446655, -0.00019259491531984653, 0.0, -90.0),
    // polar ellipsoid, inverse (main: the pole, 0.11 mm and 0.22 mm off)
    ("+proj=aeqd +lat_0=90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.999999999, 5.584744029049776e-05, -9.673060405581132e-05, 5.584694160221114e-05, -9.672974030236177e-05, 0.0, 90.0),
    ("+proj=aeqd +lat_0=-90 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", -150.0, -89.999999998, -0.00011169397415233011, -0.0001934596381311207, -0.00011169388320442228, -0.00019345948060472355, 0.0, -90.0),
    // geodesic: forward (0, 0) and inverse to the centre (main: 0.08 mm, 0.14 mm and 0.22 mm off)
    ("+proj=aeqd +lat_0=45 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 1e-09, 45.0, 7.884683509462603e-05, 0.0, 0.0, 0.0, 0.0, 45.0),
    ("+proj=aeqd +lat_0=45 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", 1e-09, 45.000000001, 7.884683509462603e-05, 0.00011113157708699083, 0.0, 0.0, 0.0, 45.0),
    ("+proj=aeqd +lat_0=-30 +ellps=GRS80 +x_0=0 +y_0=0 +lon_0=0", -2e-09, -30.000000001, -0.00019297256050258458, -0.00011085294781510523, 0.0, 0.0, 0.0, -29.999999999999996),
];

/// Checks `cases` against their exact values (forward and inverse, within `tol` metres) and
/// against the PROJ 9.8.1 outputs (see `check_proj()`). The figure comes from the definition:
/// `+R=<radius>`, or GRS80.
fn check_cases(cases: &[Case], tol: f64) -> Vec<String> {
    let mut failures = Vec::new();
    for &(def, lon, lat, x, y, px, py, plon, plat) in cases {
        let r = def
            .split_whitespace()
            .find_map(|p| p.strip_prefix("+R="))
            .map(|v| v.parse::<f64>().unwrap());
        let (ll, r) = match r {
            Some(r) => (format!("+proj=latlong +R={r}"), r),
            None => ("+proj=latlong +ellps=GRS80".to_string(), A),
        };
        let dist = |l1: f64, p1: f64, l2: f64, p2: f64| distance(l1, p1, l2, p2) * r / R;
        let ll = Proj::from_proj_string(&ll).unwrap();
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (lon.to_radians(), lat.to_radians(), 0.);
        match transform(&ll, &pj, &mut p) {
            Ok(()) => {
                let fwd = (p.0 - x).hypot(p.1 - y);
                if fwd.is_nan() || fwd >= tol {
                    failures.push(format!(
                        "{def} forward ({lon}, {lat}): ({}, {}) is {fwd:e} m off",
                        p.0, p.1
                    ));
                }
                check_proj(
                    &mut failures,
                    format!("{def} forward ({lon}, {lat})"),
                    (p.0 - px).hypot(p.1 - py),
                    fwd,
                    (px - x).hypot(py - y),
                );
            }
            Err(e) => failures.push(format!("{def} forward ({lon}, {lat}) failed: {e:?}")),
        }
        let mut q = (x, y, 0.);
        if let Err(e) = transform(&pj, &ll, &mut q) {
            failures.push(format!("{def} inverse ({x}, {y}) failed: {e:?}"));
            continue;
        }
        let (lon, lat, plon, plat) = (
            lon.to_radians(),
            lat.to_radians(),
            plon.to_radians(),
            plat.to_radians(),
        );
        let inv = dist(q.0, q.1, lon, lat);
        if inv.is_nan() || inv >= tol {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {inv:e} m off",
                q.0.to_degrees(),
                q.1.to_degrees()
            ));
        }
        check_proj(
            &mut failures,
            format!("{def} inverse ({x}, {y})"),
            dist(q.0, q.1, plon, plat),
            inv,
            dist(plon, plat, lon, lat),
        );
    }
    failures
}

#[test]
fn aeqd_near_centre() {
    // 1e-8 m: in the polar aspects the forward takes pi/2 + phi or Mp - M(phi), whose rounding
    // is about 1e-9 m; one unit in the last place of the latitude is 1.4e-9 m near the poles
    let failures = check_cases(AEQD_NEAR_CENTRE, 1e-8);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// tmerc on the Moon sphere, close to the equator, at longitudes beyond 90 degrees from the
/// central meridian (which proj4rs rejected) and close to +/-90 degrees, from OSGeo/PROJ#4945.
/// Exact values are R * atanh(cos(phi) sin(lam)) and R * atan2(sin(phi), cos(phi) cos(lam)) from a
/// 70-digit computation on the same (double radians) input; they match the PR's to 1e-7 m except
/// close to +/-90 degrees, where one unit in the last place of the longitude moves the easting by
/// about 0.01 mm. PROJ 9.8.1 does not include #4945: close to the equator its northing is 0 or
/// up to 9 mm off, and close to +/-90 degrees it is up to 0.24 m off. See `Case` for the fields.
#[rustfmt::skip]
const TMERC_MOON: &[Case] = &[
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 0.0, -1e-09, 0.0, -3.0323350424149487e-05, 0.0, 0.0, 0.0, -0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 0.0, 1e-09, 0.0, 3.0323350424149487e-05, 0.0, 0.0, 0.0, 0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 3.5, -1e-09, 106197.79433882337, -3.038001527658708e-05, 106197.79433882322, 0.0, 3.499999999999996, -0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 2.0, -3e-07, 60659.020628991224, -0.009102550154856509, 60659.02062899121, 0.0, 2.0000000000000058, -0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 1.0, -1e-06, 30324.890044860662, -0.030327969520804316, 30324.890044860545, -0.025889277458190918, 0.9999999999999967, -8.536436123899243e-07),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 0.5, 4e-06, 15161.867653671146, 0.12129802035353507, 15161.867653671059, 0.12143147500034017, 0.5000000000000061, 4.0044008845788636e-06),
    // beyond 90 degrees from the central meridian (rejected before)
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 120.0, -1e-09, 2288082.6501171775, -5458203.07628626, 2288082.650117177, -5458203.076346907, 119.99999999999999, -0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 120.0, 1e-09, 2288082.6501171775, 5458203.07628626, 2288082.650117177, 5458203.076346907, 119.99999999999999, 0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", -150.0, 1e-07, -954364.4951659868, 5458203.072845468, -954364.4951659869, 5458203.076346907, -150.0, 0.0),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", -150.0, -1e-06, -954364.4951659866, -5458203.041332518, -954364.4951659866, -5458203.050457629, -150.0, -1.0456548945759365e-06),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 179.0, -1e-09, 30324.890044860535, -5458203.076316579, 30324.890044860545, -5458203.076346907, 179.0, -0.0),
    // the acos() argument rounded to -1 and gave y = 0 instead of +/- pi * R in PROJ
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 172.622541429877, -9.247995698785671e-07, 224329.99619768397, -5458203.0480698, 224329.99619768382, 0.0, 172.622541429877, -8.467059002667616e-07),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 172.622541429877, 9.247995698785671e-07, 224329.99619768397, 5458203.0480698, 224329.99619768382, 0.0, 172.622541429877, 8.467059002667616e-07),
    // close to +/-90 degrees at the equator, where log((1 + b) / (1 - b)) lost up to 0.2 m
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", -89.998988, -0.001034, -19597452.488251183, -1383231.76970209, -19597452.707324892, -1383231.5554250323, -89.998988, -0.0010339999999999987),
    ("+proj=tmerc +R=1737400 +x_0=0 +y_0=0 +lon_0=0", 90.001468, 0.001013, 19233888.928566705, 4408778.908582411, 19233889.09228675, 4408779.145534157, 90.001468, 0.0010129999999999994),
];

#[test]
fn tmerc_sphere_moon() {
    let failures = check_cases(TMERC_MOON, 1e-6);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
