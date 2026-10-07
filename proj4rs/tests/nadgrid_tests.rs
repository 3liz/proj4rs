//! Local nadgrid tests
//!
//! $PROJ_NAGRIDS must point to a valid nadgrid repository
//!
#[cfg(feature = "local_tests")] 
mod local_tests {

use approx::assert_abs_diff_eq;
use proj4rs::nadgrids::{self, catalog};
use proj4rs::Proj;
use proj4rs::transform::transform;
use std::sync::Once;

use env_logger;

static INIT: Once = Once::new();

pub fn setup() {
    // Init setup
    INIT.call_once(|| {
        env_logger::init();
        catalog::set_builder(nadgrids::files::read_from_file);
    });
}

#[test]
fn test_wgs84_bng_nadgrid_conversion() {
    setup();

    let from = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();
    let to = Proj::from_proj_string(concat!(
        "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 ",
        "+ellps=airy +nadgrids=europe/OSTN15_NTv2_OSGBtoETRS.gsb",
    ))
    .unwrap();

    let mut v1 = vec![(-4.89328_f64.to_radians(), 51.66311_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    eprintln!("{:?}", v1[0]);

    assert_abs_diff_eq!(v1[0].0, 199999.973939543968, epsilon = 1.0e-6);
    assert_abs_diff_eq!(v1[0].1, 200000.366094537778, epsilon = 1.0e-6);
}


#[test]
#[cfg(feature = "local_tests")]
fn test_wgs84_bng_latlong_nadgrid() {
    setup();

    let from = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();
    let to = Proj::from_proj_string(concat!(
        "+proj=latlong ",
        "+nadgrids=europe/OSTN15_NTv2_OSGBtoETRS.gsb",
    ))
    .unwrap();

    let mut v1 = vec![(-9.0_f64.to_radians(), 49.0_f64.to_radians(), 0.0)];
    //let mut v1 = vec![(-4.89328_f64.to_radians(), 51.66311_f64.to_radians(), 0.0)];

    transform(&from, &to, v1.as_mut_slice()).unwrap();

    eprintln!("{:?}", (v1[0].0.to_degrees(), v1[0].1.to_degrees()));

    // Compare to output of proj 9
    // echo -9.0 49.0 | cct -z0 -t0 -d 12 +proj=longlat +nadgrids=OSTN15_NTv2_OSGBtoETRS.gsb
    assert_abs_diff_eq!(v1[0].0, -8.999464150263_f64.to_radians(), epsilon = 1.0e-10);
    assert_abs_diff_eq!(v1[0].1, 48.999301262247_f64.to_radians(), epsilon = 1.0e-10);
}

#[test]
#[cfg(feature = "local_tests")]
fn test_epsg27700_bad_point() {
    // From https://github.com/3liz/proj4rs/issues/37
    setup();

    const EPSG_27700: &str = concat!(
        "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 ",
        "+ellps=airy +nadgrids=europe/OSTN15_NTv2_OSGBtoETRS.gsb",
    );

    let epsg_4326 = Proj::from_proj_string("+proj=longlat +datum=WGS84").unwrap();
    let epsg_27700 = Proj::from_proj_string(EPSG_27700).unwrap();

    proj4rs::adaptors::transform_vertex_2d(
        &epsg_4326,
        &epsg_27700,
        (-0.03209530211282055, 0.8866271675445546),
    )
    .unwrap();

    proj4rs::adaptors::transform_vertex_2d(&epsg_4326, &epsg_27700, (-0.0321, 0.8866271675445546))
        .unwrap();

    proj4rs::adaptors::transform_vertex_2d(
        &epsg_4326,
        &epsg_27700,
        (-0.03209530211282055, 0.8866272),
    )
    .unwrap();
}


}
