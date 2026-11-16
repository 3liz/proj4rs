//! Parity with PROJ 9.8.1 at ordinary points, for the projections whose precision or formulas were
//! changed (see PRECISIONS.md). `fixtures/proj_parity/proj_9_8_1.csv` holds PROJ's forward and inverse
//! outputs, generated with pyproj by `fixtures/proj_parity/generate.py`. The points avoid the
//! places where PROJ itself loses precision; those are covered, against exact values, by
//! `tests/precision_tests.rs`.

use proj4rs::Proj;
use proj4rs::transform::transform;

const DATA: &str = include_str!("../fixtures/proj_parity/proj_9_8_1.csv");

/// Tolerance in metres, on the projected coordinates and on the inverse (as a distance)
const TOL: f64 = 1e-6;

/// Distance in metres between two points given in radians, on a sphere of radius r
fn distance(l1: f64, p1: f64, l2: f64, p2: f64, r: f64) -> f64 {
    let a = ((p2 - p1) / 2.).sin();
    let b = ((l2 - l1) / 2.).sin();
    2. * r * (a * a + p1.cos() * p2.cos() * b * b).sqrt().asin()
}

#[test]
fn parity_with_proj_9_8_1() {
    let mut failures = Vec::new();
    let mut count = 0;
    let mut current: Option<(String, Proj, Proj, f64)> = None;
    for line in DATA
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let v: Vec<&str> = line.split(';').collect();
        assert_eq!(v.len(), 7, "bad fixture line: {line}");
        let def = v[0];
        let n: Vec<f64> = v[1..].iter().map(|s| s.parse().unwrap()).collect();
        let (lon, lat, x, y, inv_lon, inv_lat) = (n[0], n[1], n[2], n[3], n[4], n[5]);

        if current.as_ref().is_none_or(|c| c.0 != def) {
            // geographic CRS on the figure of the definition; distances on a sphere of
            // radius R, or a for an ellipsoid (close enough for the distances checked here)
            let figure: Vec<&str> = def
                .split_whitespace()
                .filter(|p| p.starts_with("+R=") || p.starts_with("+ellps="))
                .collect();
            let ll = format!("+proj=latlong {}", figure.join(" "));
            let r = def
                .split_whitespace()
                .find_map(|p| p.strip_prefix("+R="))
                .map_or(6378137., |v| v.parse().unwrap());
            current = Some((
                def.to_string(),
                Proj::from_proj_string(&ll).unwrap(),
                Proj::from_proj_string(def).unwrap(),
                r,
            ));
        }
        let (_, ll, pj, r) = current.as_ref().unwrap();
        count += 1;

        let mut p = (lon.to_radians(), lat.to_radians(), 0.);
        match transform(ll, pj, &mut p) {
            Ok(()) => {
                let d = (p.0 - x).hypot(p.1 - y);
                if d.is_nan() || d >= TOL {
                    failures.push(format!(
                        "{def} forward ({lon}, {lat}): ({}, {}), PROJ ({x}, {y}): {d:e} m",
                        p.0, p.1
                    ));
                }
            }
            Err(e) => failures.push(format!("{def} forward ({lon}, {lat}) failed: {e:?}")),
        }

        let mut q = (x, y, 0.);
        match transform(pj, ll, &mut q) {
            Ok(()) => {
                let d = distance(q.0, q.1, inv_lon.to_radians(), inv_lat.to_radians(), *r);
                if d.is_nan() || d >= TOL {
                    failures.push(format!(
                        "{def} inverse ({x}, {y}): ({}, {}), PROJ ({inv_lon}, {inv_lat}): {d:e} m",
                        q.0.to_degrees(),
                        q.1.to_degrees()
                    ));
                }
            }
            Err(e) => failures.push(format!("{def} inverse ({x}, {y}) failed: {e:?}")),
        }
    }
    assert!(count > 0, "empty fixture");
    assert!(
        failures.is_empty(),
        "{} of {count} points differ from PROJ 9.8.1:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
