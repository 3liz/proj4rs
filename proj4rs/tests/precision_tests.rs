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
    ("+proj=stere +lat_0=90 +R=6371000 +k_0=0.994 +x_0=0 +y_0=0 +lon_0=0", 30.0, 89.9999999, 0.0055263875261430016, -0.0095719839775945559),
    // ortho sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=ortho +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 10.0, 89.9999999, 0.0019308795231307598, 4504977.2951962777),
    // aeqd sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=aeqd +lat_0=-30 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -120.0, -89.9999999, -0.011644304801302344, -6671695.6042332702),
    // laea sphere, oblique: latitude from asin() lost precision near the pole
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 60.0, 89.9999999, 0.01042317987605592, 4876152.2900594383),
    // laea sphere, equatorial: points within 1e-10 * R of the centre came back on the equator
    ("+proj=laea +lat_0=0 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 1e-09, -2e-09, 0.00011119492664455874, -0.00022238985328911749),
    // gnom sphere, polar: latitude from asin() lost precision near the pole
    ("+proj=gnom +lat_0=-90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 45.0, -89.9999999, 0.0078626681995985214, 0.0078626681995985214),
    // gnom: the centre itself came back with an undefined latitude
    ("+proj=gnom +lat_0=90 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 90.0, 0.0, 0.0),
    // tpers: latitude from asin() lost precision near the pole
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.9999999, 0.0067512397377098711, 2777408.3169169292),
    // tpers: the centre itself came back with latitude 0 instead of lat_0
    ("+proj=tpers +lat_0=45 +R=6371000 +h=3000000 +x_0=0 +y_0=0 +lon_0=0", 0.0, 45.0, 0.0, 0.0),
    // moll: near the poles the Newton iteration stopped early and the inverse clamped sin(theta)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", -179.0, -89.99999, -591.0783918084689, -9009954.6009776294),
    // moll: near the poles the Newton iteration stopped early and the inverse clamped sin(theta)
    ("+proj=moll +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 100.0, 89.99, 33021.103249211675, 9009905.5922611956),
    // tmerc sphere: close to the equator acos() lost half the digits
    ("+proj=tmerc +R=6371000 +k=0.999 +x_0=0 +y_0=0 +lon_0=0 +lat_0=0", 2.0, 1e-07, 222212.59457055129, 0.011115144217327468),
    // laea sphere, oblique: the inverse gave longitude 0 instead of 90
    ("+proj=laea +lat_0=45 +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 90.0, -10.0, 9473745.8826114854, -1181205.6198054144),
    // eqearth: the inverse gave NaN at the poles
    ("+proj=eqearth +R=6371000 +x_0=0 +y_0=0 +lon_0=0", 30.0, 90.0, 1702743.8797029482, 8392918.138591881),
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
        if !(fwd < 1e-6) {
            failures.push(format!("{def} forward ({lon}, {lat}): ({}, {}) is {fwd:e} m off", p.0, p.1));
        }
        let mut q = (x, y, 0.);
        transform(&pj, &ll, &mut q).unwrap();
        let inv = distance(q.0, q.1, lon.to_radians(), lat.to_radians());
        if !(inv < 1e-6) {
            failures.push(format!("{def} inverse ({x}, {y}): ({}, {}) is {inv:e} m off", q.0.to_degrees(), q.1.to_degrees()));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

