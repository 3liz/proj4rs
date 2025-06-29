//!
//! Tests from proj4js
//!
//! Note: projection results may differs from proj by 10^-4 due to difference
//! in math functions implementations (asinh, log1py...)
//!
use approx::assert_abs_diff_eq;
use proj4rs::{proj, transform};

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

    let from = proj::Proj::from_user_string(sweref99tm).unwrap();
    let to = proj::Proj::from_user_string(rt90).unwrap();

    let mut inp = (319180., 6399862., 0.);

    transform::transform(&from, &to, &mut inp).unwrap();
    assert_abs_diff_eq!(inp.0, 1271137.92755580, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 6404230.29136189, epsilon = 1.0e-6);
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

    let from = proj::Proj::from_user_string(epsg2154).unwrap();
    let to = proj::Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (489353.59, 6587552.2, 0.);
    transform::transform(&from, &to, &mut inp).unwrap();
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

    let projection = proj::Proj::from_user_string(wgs84);
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

    let from = proj::Proj::from_user_string(epsg3044).unwrap();
    let to = proj::Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (580900., 5625000., 0.);
    transform::transform(&from, &to, &mut inp).unwrap();
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

    let from = proj::Proj::from_user_string(epsg3044).unwrap();
    let to = proj::Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (580900., 5625000., 0.);
    transform::transform(&from, &to, &mut inp).unwrap();
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

    let from = proj::Proj::from_user_string(epsg3844).unwrap();
    let to = proj::Proj::from_user_string(epsg3857).unwrap();

    let mut inp = (505000., 500000., 0.);
    transform::transform(&from, &to, &mut inp).unwrap();
    // Compare results from cs2cs output
    assert_abs_diff_eq!(inp.0, 2790174.2500622645, epsilon = 1.0e-6);
    assert_abs_diff_eq!(inp.1, 5780346.2980352566, epsilon = 1.0e-6);
}
