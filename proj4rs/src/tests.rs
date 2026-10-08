//!
//! Unit tests
//!
use env_logger;
use std::sync::Once;

static INIT: Once = Once::new();

pub fn setup() {
    // Init setup
    INIT.call_once(|| {
        env_logger::init();
    });
}

pub(crate) mod utils {
    use crate::math::adjlon;
    use crate::proj::{Proj, ProjData};
    use approx::assert_abs_diff_eq;

    pub(crate) fn scale(d: &ProjData, xyz: (f64, f64, f64)) -> (f64, f64, f64) {
        (xyz.0 * d.ellps.a + d.x0, xyz.1 * d.ellps.a + d.y0, xyz.2)
    }

    pub(crate) fn descale(d: &ProjData, xyz: (f64, f64, f64)) -> (f64, f64, f64) {
        (
            (xyz.0 - d.x0) * d.ellps.ra,
            (xyz.1 - d.y0) * d.ellps.ra,
            xyz.2,
        )
    }

    pub(crate) fn to_deg(lam: f64, phi: f64, z: f64) -> (f64, f64, f64) {
        (lam.to_degrees(), phi.to_degrees(), z)
    }

    pub(crate) fn to_rad(lpz: (f64, f64, f64)) -> (f64, f64, f64) {
        (lpz.0.to_radians(), lpz.1.to_radians(), lpz.2)
    }

    pub(crate) fn test_proj_forward(
        p: &Proj,
        inputs: &[((f64, f64, f64), (f64, f64, f64))],
        prec: f64,
    ) {
        let d = p.data();
        inputs.iter().for_each(|(input, expect)| {
            let (lam, phi, z) = to_rad(*input);
            let out = scale(
                d,
                p.projection()
                    .forward(adjlon(lam - d.lam0), phi, z)
                    .unwrap(),
            );
            println!("{:?}", out);
            assert_abs_diff_eq!(out.0, expect.0, epsilon = prec);
            assert_abs_diff_eq!(out.1, expect.1, epsilon = prec);
            assert_abs_diff_eq!(out.2, expect.2, epsilon = prec);
        })
    }

    pub(crate) fn test_proj_inverse(
        p: &Proj,
        inputs: &[((f64, f64, f64), (f64, f64, f64))],
        prec: f64,
    ) {
        let d = p.data();
        inputs.iter().for_each(|(expect, input)| {
            let (x, y, z) = descale(d, *input);
            let (lam, phi, z) = p.projection().inverse(x, y, z).unwrap();
            let out = to_deg(adjlon(lam + d.lam0), phi, z);
            println!("{:?}", out);
            assert_abs_diff_eq!(out.0, expect.0, epsilon = prec);
            assert_abs_diff_eq!(out.1, expect.1, epsilon = prec);
            assert_abs_diff_eq!(out.2, expect.2, epsilon = prec);
        })
    }
}

use crate::proj::Proj;
use crate::transform::transform;
use approx::assert_abs_diff_eq;


#[test]
fn test_transform_array() {
    let mut data: Vec<(f64, f64, f64)> = (1..=1_000)
        .map(|_| (2.0f64.to_radians(), 1.0f64.to_radians(), 0.0f64))
        .collect();

    let from = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
    let to = Proj::from_proj_string("+proj=etmerc +ellps=GRS80").unwrap();

    transform(&from, &to, data.as_mut_slice()).unwrap();

    // Check values
    data.iter().for_each(|(x, y, _)| {
        assert_abs_diff_eq!(*x, 222650.79679758527, epsilon = 1.0e-10);
        assert_abs_diff_eq!(*y, 110642.22941193319, epsilon = 1.0e-10);
    });
}

#[test]
fn test_utm33_grs80() {
    let from = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
    let to = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=33").unwrap();

    let mut v1 = vec![(
        13.393921852111816_f64.to_radians(),
        52.5200080871582_f64.to_radians(),
        0.0,
    )];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v1[0].0, 391027.67777461524, epsilon = 1.0e-10);
    assert_abs_diff_eq!(v1[0].1, 5820089.724404063, epsilon = 1.0e-10);
}

#[test]
fn test_wgs84_bng_conversion() {
    //crate::nadgrids::catalog::files::

    let from = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();
    let to = Proj::from_proj_string(concat!(
        "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 ",
        "+ellps=airy ", //+nadgrids=OSTN15_NTv2_OSGBtoETRS.gsb",
    ))
    .unwrap();

    let mut v1 = vec![(-4.89328_f64.to_radians(), 51.66311_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v1[0].0, 199925.978901151626, epsilon = 1.0e-8);
    assert_abs_diff_eq!(v1[0].1, 200052.051949012151, epsilon = 1.0e-8);
}

#[test]
fn test_prime_meridian_parameter_by_name() {
    let from = Proj::from_proj_string("+proj=merc +ellps=WGS84 +pm=paris").unwrap();
    let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

    let mut v1 = vec![(0.0_f64.to_radians(), 0.0_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v1[0].0.to_degrees(), 2.337, epsilon = 1.0e-3);
}

#[test]
fn test_prime_meridian_parameter_dms() {
    let from = Proj::from_proj_string("+proj=merc +ellps=WGS84 +pm=2d20'14.025\"E").unwrap();
    let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

    let mut v1 = vec![(0.0_f64.to_radians(), 0.0_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v1[0].0.to_degrees(), 2.337, epsilon = 1.0e-3);
}

#[test]
fn test_prime_meridian_for_latlong_geocent() {
    let from = Proj::from_proj_string("+proj=latlong +ellps=WGS84 +pm=paris").unwrap();
    let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

    let mut v1 = vec![(0.0_f64.to_radians(), 0.0_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v1[0].0.to_degrees(), 2.337, epsilon = 1.0e-3);
}

#[test]
// Regression: geocentric output must be converted from meters to
// the target units (divided by `to_meter`), not multiplied.
fn test_geocent_units_scaling() {
    let from = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();
    let to = Proj::from_proj_string("+proj=geocent +datum=WGS84 +units=km").unwrap();

    let mut v = vec![
        (0.0_f64.to_radians(), 0.0_f64.to_radians(), 0.0),
        (0.0_f64.to_radians(), 90.0_f64.to_radians(), 0.0),
    ];

    transform(&from, &to, v.as_mut_slice()).unwrap();

    // Equator / Greenwich: X = a
    assert_abs_diff_eq!(v[0].0, 6378.137, epsilon = 1.0e-9);
    assert_abs_diff_eq!(v[0].1, 0., epsilon = 1.0e-9);
    assert_abs_diff_eq!(v[0].2, 0., epsilon = 1.0e-9);
    // North pole: Z = b
    assert_abs_diff_eq!(v[1].0, 0., epsilon = 1.0e-9);
    assert_abs_diff_eq!(v[1].1, 0., epsilon = 1.0e-9);
    assert_abs_diff_eq!(v[1].2, 6356.752314245, epsilon = 1.0e-9);

    // Round trip back to geographic
    transform(&to, &from, v.as_mut_slice()).unwrap();

    assert_abs_diff_eq!(v[0].0, 0., epsilon = 1.0e-12);
    assert_abs_diff_eq!(v[0].1, 0., epsilon = 1.0e-12);
    assert_abs_diff_eq!(v[0].2, 0., epsilon = 1.0e-6);
    assert_abs_diff_eq!(v[1].1, 90.0_f64.to_radians(), epsilon = 1.0e-12);
    assert_abs_diff_eq!(v[1].2, 0., epsilon = 1.0e-6);
}
