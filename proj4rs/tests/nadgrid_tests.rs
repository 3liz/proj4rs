//! Local nadgrid tests for ntv2 grids
//!
//! These tests use the OSGeo datum repository: https://github.com/OSGeo/proj-datumgrid
//!
//! Clone that repo and set the `PROJ_NAGRIDS` environment variable to the path 
//! of the repo.
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
fn test_nadgrid_subgrid_lookup() {
    // Check that the subgrid hierarchy is correctly walked
    //
    // Regression: subgrid search must not stop on a grandchild grid and skip
    // the remaining sibling subgrids.
    //
    // SK83-98.GSB hierarchy: SKcsrs5m (root) with 15 children; 'fortq30s'
    // has one child 'fqsub03s'.
    setup();

    let from = Proj::from_proj_string(
        "+proj=latlong +ellps=GRS80 +nadgrids=north-america/SK83-98.GSB",
    )
    .unwrap();
    let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

    // Input at the centre of each subgrid, expected output from proj 9.8.1:
    // echo <lon> <lat> | cct -z0 -t0 -d 12 \
    //    +proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad \
    //    +step +proj=hgridshift +grids=SK83-98.GSB \
    //    +step +proj=unitconvert +xy_in=rad +xy_out=deg
    let cases: [(&str, (f64, f64), (f64, f64)); 4] = [
        // Sibling after the grandchild
        ("chous30s", (-102.875, 53.75), (-102.875014238888, 53.750006552778)),
        // Sibling after the grandchild
        ("grass30s", (-106.70833, 49.08333), (-106.708331800553, 49.083323900542)),
        // Sibling after the grandchild
        ("batt30s", (-108.16667, 52.79167), (-108.166670230549, 52.791669594437)),
        // The grandchild itself
        ("fqsub03s", (-103.67917, 50.75417), (-103.679159748695, 50.754174386961)),
    ];

    for (name, (lon, lat), (exp_lon, exp_lat)) in cases {
        let mut v = (lon.to_radians(), lat.to_radians(), 0.0);
        transform(&from, &to, &mut v).unwrap();

        let (out_lon, out_lat) = (v.0.to_degrees(), v.1.to_degrees());
        eprintln!("{name}: ({out_lon:.12}, {out_lat:.12})");

        // 1e-8 deg ~ 1 mm
        assert_abs_diff_eq!(out_lon, exp_lon, epsilon = 1.0e-8);
        assert_abs_diff_eq!(out_lat, exp_lat, epsilon = 1.0e-8);
    }
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
