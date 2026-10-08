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
    use proj4rs::Proj;
    use proj4rs::nadgrids::{self, catalog};
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
        // Regression test: subgrid search must not stop on a grandchild grid and skip
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
            (
                "chous30s",
                (-102.875, 53.75),
                (-102.875014238888, 53.750006552778),
            ),
            // Sibling after the grandchild
            (
                "grass30s",
                (-106.70833, 49.08333),
                (-106.708331800553, 49.083323900542),
            ),
            // Sibling after the grandchild
            (
                "batt30s",
                (-108.16667, 52.79167),
                (-108.166670230549, 52.791669594437),
            ),
            // The grandchild itself
            (
                "fqsub03s",
                (-103.67917, 50.75417),
                (-103.679159748695, 50.754174386961),
            ),
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
    fn test_nadgrid_null_fallback() {
        // A trailing '@null' in a grid list acts as a world-wide
        // zero-shift fallback: points outside the other grids must pass through
        // unchanged instead of failing.
        setup();

        let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

        // ntf_r93.gsb covers lat [41, 52], lon [-5.5, 10]
        let from =
            Proj::from_proj_string("+proj=latlong +ellps=clrk80ign +nadgrids=ntf_r93.gsb,@null")
                .unwrap();

        // Expected output from proj 9.8.1:
        // cs2cs -d 12 +proj=longlat +ellps=clrk80ign +nadgrids=ntf_r93.gsb,@null \
        //    +to +proj=longlat +datum=WGS84
        let cases: [(&str, (f64, f64), (f64, f64)); 5] = [
            // Inside the grid: shifted
            ("Paris", (2.35, 48.85), (2.349295593686, 48.849933562569)),
            ("Brest", (-4.5, 48.4), (-4.500970456948, 48.399916990002)),
            // Outside the grid: '@null' fallback, unchanged
            ("Prague", (15.0, 50.0), (15.0, 50.0)),
            ("Atlantic", (-10.0, 40.0), (-10.0, 40.0)),
            ("East of grid", (10.5, 47.0), (10.5, 47.0)),
        ];

        for (name, (lon, lat), (exp_lon, exp_lat)) in cases {
            let mut v = (lon.to_radians(), lat.to_radians(), 0.0);
            transform(&from, &to, &mut v).unwrap_or_else(|e| panic!("{name}: {e:?}"));

            let (out_lon, out_lat) = (v.0.to_degrees(), v.1.to_degrees());
            eprintln!("{name}: ({out_lon:.12}, {out_lat:.12})");

            // 1e-8 deg ~ 1 mm
            assert_abs_diff_eq!(out_lon, exp_lon, epsilon = 1.0e-8);
            assert_abs_diff_eq!(out_lat, exp_lat, epsilon = 1.0e-8);
        }

        // Without '@null', points outside the grid are still an error
        let from =
            Proj::from_proj_string("+proj=latlong +ellps=clrk80ign +nadgrids=ntf_r93.gsb").unwrap();

        let mut v = (15.0_f64.to_radians(), 50.0_f64.to_radians(), 0.0);
        assert!(transform(&from, &to, &mut v).is_err());
    }

    #[test]
    #[cfg(feature = "local_tests")]
    fn test_nadgrid_identical_datums() {
        // Two CRS with the same nadgrids list have identical datums,
        // so no grid shift must be applied between them. A point outside the grid
        // must then pass through unchanged instead of failing.
        setup();

        fn latlong(nadgrids: &str) -> Proj {
            Proj::from_proj_string(&format!(
                "+proj=latlong +ellps=clrk80ign +nadgrids={nadgrids}"
            ))
            .unwrap()
        }

        // ntf_r93.gsb covers lat [41, 52], lon [-5.5, 10]
        let inside = (2.35_f64, 48.85_f64); // Paris
        let outside = (15.0_f64, 50.0_f64); // Prague

        // Results checked against proj 9.8.1:
        // cs2cs -d 15 +proj=longlat +ellps=clrk80ign +nadgrids=<src> \
        //    +to +proj=longlat +ellps=clrk80ign +nadgrids=<dst>
        let identical = [
            ("ntf_r93.gsb", "ntf_r93.gsb"),
            ("ntf_r93.gsb,@null", "ntf_r93.gsb,@null"),
        ];

        for (src, dst) in identical {
            let (from, to) = (latlong(src), latlong(dst));
            for (lon, lat) in [inside, outside] {
                let mut v = (lon.to_radians(), lat.to_radians(), 0.0);
                transform(&from, &to, &mut v)
                    .unwrap_or_else(|e| panic!("{src} -> {dst} ({lon}, {lat}): {e:?}"));

                // No shift, not even a forward/inverse round trip
                assert_abs_diff_eq!(v.0.to_degrees(), lon, epsilon = 1.0e-12);
                assert_abs_diff_eq!(v.1.to_degrees(), lat, epsilon = 1.0e-12);
            }
        }

        // Different grid lists are not identical datums: the source grid is
        // applied and fails outside of its area
        let (from, to) = (latlong("ntf_r93.gsb"), latlong("ntf_r93.gsb,@null"));
        let mut v = (outside.0.to_radians(), outside.1.to_radians(), 0.0);
        assert!(transform(&from, &to, &mut v).is_err());
    }

    #[test]
    #[cfg(feature = "local_tests")]
    fn test_nadgrid_multiple_roots() {
        // All top-level grids of a file must be used, not only
        // the first one.
        //
        // ntv2_0.gsb (NAD27 -> NAD83, Canada) has 4 root grids, in file order:
        // CAeast, CAwest, CAnorth and CAarctic.
        setup();

        let from = Proj::from_proj_string(
            "+proj=latlong +ellps=clrk66 +nadgrids=north-america/ntv2_0.gsb",
        )
        .unwrap();
        let to = Proj::from_proj_string("+proj=latlong +datum=WGS84").unwrap();

        // Expected output from proj 9.8.1:
        // cs2cs -d 12 +proj=longlat +ellps=clrk66 +nadgrids=ntv2_0.gsb \
        //    +to +proj=longlat +datum=WGS84
        let cases: [(&str, (f64, f64), (f64, f64)); 7] = [
            // First root
            (
                "Montreal (CAeast)",
                (-73.57, 45.50),
                (-73.569587678433, 45.500039887888),
            ),
            // Other roots: were failing with PointOutsideNadShiftArea
            (
                "Calgary (CAwest)",
                (-114.07, 51.05),
                (-114.070992067239, 51.050056225001),
            ),
            (
                "Vancouver (CAwest)",
                (-123.12, 49.28),
                (-123.121321768916, 49.279828606221),
            ),
            (
                "Yellowknife (CAnorth)",
                (-114.37, 62.45),
                (-114.371273380006, 62.450197115009),
            ),
            (
                "Iqaluit (CAnorth)",
                (-68.52, 63.75),
                (-68.518895049399, 63.750296782788),
            ),
            (
                "Eureka (CAarctic)",
                (-85.93, 79.99),
                (-85.929649486567, 79.991017875816),
            ),
            (
                "Alert (CAarctic)",
                (-62.35, 82.50),
                (-62.345059231252, 82.500971743320),
            ),
        ];

        for (name, (lon, lat), (exp_lon, exp_lat)) in cases {
            let mut v = (lon.to_radians(), lat.to_radians(), 0.0);
            transform(&from, &to, &mut v).unwrap_or_else(|e| panic!("{name}: {e:?}"));

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

        proj4rs::adaptors::transform_vertex_2d(
            &epsg_4326,
            &epsg_27700,
            (-0.0321, 0.8866271675445546),
        )
        .unwrap();

        proj4rs::adaptors::transform_vertex_2d(
            &epsg_4326,
            &epsg_27700,
            (-0.03209530211282055, 0.8866272),
        )
        .unwrap();
    }
}
