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

#[test]
// Regression test for issue #34: a standalone `+nadgrids=@null` means a zero
// datum shift (identity to WGS84), NOT an unknown datum: the +towgs84
// of the other side must still be applied.
fn test_nadgrids_null_keeps_other_side_towgs84() {
    let from = Proj::from_proj_string(concat!(
        "+proj=sterea +lat_0=52.1561605555556 +lon_0=5.38763888888889 ",
        "+k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel ",
        "+towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 ",
        "+units=m +no_defs",
    ))
    .unwrap();
    let to = Proj::from_proj_string(concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 ",
        "+k=1 +units=m +nadgrids=@null +no_defs",
    ))
    .unwrap();

    // Amersfoort tower, RD coordinates (155000, 463000).
    let mut v = vec![(155000.0_f64, 463000.0_f64, 0.0)];
    transform(&from, &to, v.as_mut_slice()).unwrap();

    // Check against cs2cs (proj 9.4) output
    assert_abs_diff_eq!(v[0].0, 599700.751406375, epsilon = 1.0e-8);
    assert_abs_diff_eq!(v[0].1, 6828231.372408830, epsilon = 1.0e-8);
}

#[test]
fn test_transform_with_datum() {
    //EPSG:3006 Definition - Sweden coordinate reference system
    let sweref99tm = concat!(
        "+proj=utm +zone=33 +ellps=GRS80 ",
        "+towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    );
    // EPSG:3021 Definition - Sweden coordinate reference system
    let rt90 = concat!(
        "+proj=tmerc +lon_0=15.808277777799999 +lat_0=0.0 +k=1.0 ",
        "+x_0=1500000.0 +y_0=0.0 +ellps=bessel ",
        "+units=m +towgs84=414.1,41.3,603.1,-0.855,2.141,-7.023,0 ",
        "+no_defs"
    );

    let from = Proj::from_user_string(sweref99tm).unwrap();
    let to = Proj::from_user_string(rt90).unwrap();

    let mut inp = (319180., 6399862., 0.);

    transform(&from, &to, &mut inp).unwrap();
    // Check against cs2cs (proj 9.4) output: proj4js gives
    // (1271137.92755580, 6404230.29136189) because it does not handle
    // +towgs84=0,0,0,0,0,0,0 as a 3 parameters transformation
    assert_abs_diff_eq!(inp.0, 1271137.927561178, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 6404230.291456630, epsilon = 1.0e-6);
}


#[test]
fn test_transform_null_datum() {
    // Test when nadgrid list is empty
    // ESPG:2154 definition
    let epsg2154 = concat!(
        "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 ",
        "+x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 ",
        "+units=m +no_defs +type=crs"
    );
    // ESPG:3857 definition
    let epsg3857 = concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 ",
        "+units=m +nadgrids=@null +wktext +no_defs +type=crs",
    );

    let from = Proj::from_user_string(epsg2154).unwrap();
    let to = Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (489353.59, 6587552.2, 0.);
    transform(&from, &to, &mut inp).unwrap();
    // Check against cs2cs output
    assert_abs_diff_eq!(inp.0, 28943.07106251, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 5837421.86634143, epsilon = 1.0e-6);
}

#[test]
fn test_longlat_alias() {
    let wgs84 = concat!(
        "+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 ",
        "+datum=WGS84 +units=degrees",
    );

    let projection = Proj::from_user_string(wgs84);
    assert!(projection.is_ok());
}

#[test]
fn test_transform_epsg3044() {
    // ESPG:3044 definition
    let epsg3044 = concat!("+proj=utm +zone=32 +ellps=GRS80 +units=m  +towgs84=0,0,0,0,0,0,0 ",);
    // ESPG:3857 definition
    let epsg3857 = concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 ",
        "+units=m +nadgrids=@null",
    );

    let from = Proj::from_user_string(epsg3044).unwrap();
    let to = Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (580900., 5625000., 0.);
    transform(&from, &to, &mut inp).unwrap();
    assert_abs_diff_eq!(inp.0, 1129592.3568078864, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 6580906.077194334, epsilon = 1.0e-6);
}

#[test]
fn test_axis_denormalize() {
    // ESPG:3044 definition
    let epsg3044 = concat!("+proj=utm +zone=32 +ellps=GRS80 +units=m  +towgs84=0,0,0,0,0,0,0 ",);
    // ESPG:3857 definition
    let epsg3857 = concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 ",
        "+units=m +nadgrids=@null +axis=neu",
    );

    let from = Proj::from_user_string(epsg3044).unwrap();
    let to = Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (580900., 5625000., 0.);
    transform(&from, &to, &mut inp).unwrap();
    assert_abs_diff_eq!(inp.0, 6580906.077194334, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 1129592.3568078864, epsilon = 1.0e-6);
}

#[test]
fn test_transform_epsg3844() {
    // ESPG:3844 definition
    let epsg3844 = concat!(
        "+proj=sterea +lat_0=46 +lon_0=25 +k=0.99975 +x_0=500000 +y_0=500000 ",
        "+ellps=krass ",
        //"+towgs84=2.329,-147.042,-92.08,0.309,-0.325,-0.497,5.69 ",
        //"+towgs84=44.107,-116.147,-54.648 ",
        //"+towgs84=28,-121,-77 ",
        "+units=m +no_defs +type=crs"
    );
    // ESPG:3857 definition
    let epsg3857 = concat!(
        "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 ",
        "+units=m",
    );

    // ESPG:3857 definition 2
    //let epsg3857 = concat!(
    //   "+proj=webmerc +ellps=WGS84 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 ",
    //   "+units=m +towgs84=0,0,0",
    //);

    let from = Proj::from_user_string(epsg3844).unwrap();
    let to = Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (505000., 500000., 0.);
    transform(&from, &to, &mut inp).unwrap();
    // Compare results from cs2cs output
    assert_abs_diff_eq!(inp.0, 2790174.2500622645, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 5780346.2980352566, epsilon = 1.0e-6);
}

#[test]
fn test_clark_1866() {
    // Test for https://github.com/3liz/proj4rs/issues/47 regression
    let clark =
        concat!("+proj=tmerc +lat_0=40 +lon_0=-74 +ellps=clrk66 +towgs84=-8,160,176 +units=m",);

    let from = Proj::from_user_string(clark).unwrap();
    let to = Proj::from_user_string("+proj=longlat +datum=WGS84").unwrap();

    let mut inp = (10000.0, 20000.0, 30.0);
    transform(&from, &to, &mut inp).unwrap();

    // Compare results from cs2cs output
    assert_abs_diff_eq!(inp.0.to_degrees(), -73.88215840173, epsilon = 1.0e-8);
    assert_abs_diff_eq!(inp.1.to_degrees(), 40.18007067501, epsilon = 1.0e-8);
}
