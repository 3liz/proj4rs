//! Precision of the spherical projections close to the poles, the equator and the projection
//! centres. Expected values are from a 50-digit computation on the same input, on a sphere of
//! radius 6371000 m. The inverse is checked as a great-circle distance in metres, because close
//! to the poles the longitude is ill-conditioned.

use proj4rs::Proj;
use proj4rs::transform::transform;

const R: f64 = 6371000.;

#[rustfmt::skip]
const CASES: &[(&str, f64, f64, f64, f64)] = &[
    // stere sphere, polar: latitude from asin() lost precision near the pole
    ("+proj=stere +lat_0=90 +R=6371000 +k_0=0.994 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.9999999, 0.005526387526143001, -0.009571983977594557),
    // ortho sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=ortho +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.9999999, 0.0019308795231307598, 4504977.295196278),
    // aeqd sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=aeqd +lat_0=-30 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -120.0, -89.9999999, -0.011644304801302343, -6671695.6042332705),
    // laea sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.9999999, 0.01042317987605592, 4876152.290059438),
    // laea sphere, equatorial: points within 1e-10 * R of the centre came back on the equator
    ("+proj=laea +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 1e-09, -2e-09, 0.00011119492664455874, -0.00022238985328911748),
    // gnom sphere, polar: latitude from asin() lost precision near the pole
    ("+proj=gnom +lat_0=-90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 45.0, -89.9999999, 0.007862668199598521, 0.007862668199598521),
    // gnom: the centre itself came back with an undefined latitude
    ("+proj=gnom +lat_0=90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 90.0, 0.0, 0.0),
    // tpers: latitude from asin() lost precision near the pole
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.9999999, 0.0067512397377098715, 2777408.316916929),
    // tpers: the centre itself came back with latitude 0 instead of lat_0
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 45.0, 0.0, 0.0),
    // moll: near the poles the Newton iteration stopped early and the inverse clamped sin(theta)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -179.0, -89.99999, -591.078391808469, -9009954.60097763),
    // moll: near the poles the Newton iteration stopped early and the inverse clamped sin(theta)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.99, 33021.10324921167, 9009905.592261195),
    // tmerc sphere: close to the equator acos() lost half the digits
    ("+proj=tmerc +R=6371000 +k=0.999 +x_0=0 +y_0=0 +lon_0=0 +lat_0=0", 2.0, 1e-07, 222212.5945705513, 0.011115144217327467),
    // laea sphere, oblique: the inverse gave longitude 0 instead of 90
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 90.0, -10.0, 9473745.882611485, -1181205.6198054145),
    // eqearth: the inverse gave NaN at the poles
    ("+proj=eqearth +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 30.0, 90.0, 1702743.8797029483, 8392918.138591882),
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
    for &(def, lon, lat, x, y) in CASES {
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
        let mut q = (x, y, 0.);
        transform(&pj, &ll, &mut q).unwrap();
        let inv = distance(q.0, q.1, lon.to_radians(), lat.to_radians());
        if inv.is_nan() || inv >= 1e-6 {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {inv:e} m off",
                q.0.to_degrees(),
                q.1.to_degrees()
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// gnom in equatorial and oblique mode, close to 90 degrees from the centre, where the projected
/// coordinates are very large and cos(z) computed as sqrt(1 - sin(z)^2) lost its digits: (x, y) is
/// the exact forward image rounded to double, and the expected (lon, lat) is the exact inverse of
/// those doubles, from a 50-digit computation.
#[rustfmt::skip]
const GNOM_FAR: &[(&str, f64, f64, f64, f64)] = &[
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -0.011119492664455874, -3650314329490979.0, -1e-07, -89.9999999),  // from (-1e-07, -89.9999999)
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 1123379.1940936304, 3706626311913009.5, 10.000000000000002, 89.9999999),  // from (10, 89.9999999)
    ("+proj=gnom +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 36503141116197.78, 0.0, 89.99999),  // from (0, 89.99999)
    ("+proj=gnom +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 516232372521323.0, 6371000.0, 90.0, 1e-06),  // from (90, 1e-06)
    ("+proj=gnom +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, -3650314070120549.0, 0.0, -44.9999999),  // from (0, -44.9999999)
    ("+proj=gnom +lat_0=-30 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 3650314070120549.0, 0.0, 59.9999999),  // from (0, 59.9999999)
];

#[test]
fn gnom_inverse_far_from_the_centre() {
    let ll = Proj::from_proj_string("+proj=latlong +R=6371000").unwrap();
    let mut failures = Vec::new();
    for &(def, x, y, lon, lat) in GNOM_FAR {
        let pj = Proj::from_proj_string(def).unwrap();
        let mut p = (x, y, 0.);
        transform(&pj, &ll, &mut p).unwrap();
        let d = distance(p.0, p.1, lon.to_radians(), lat.to_radians());
        if d.is_nan() || d >= 1e-6 {
            failures.push(format!(
                "{def} inverse ({x}, {y}): ({}, {}) is {d:e} m off",
                p.0.to_degrees(),
                p.1.to_degrees()
            ));
        }
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
