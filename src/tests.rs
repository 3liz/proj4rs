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
    use crate::errors::Result;
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
            let out = scale(d, p.projection().forward(lam - d.lam0, phi, z).unwrap());
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
            let out = to_deg(lam + d.lam0, phi, z);
            assert_abs_diff_eq!(out.0, expect.0, epsilon = prec);
            assert_abs_diff_eq!(out.1, expect.1, epsilon = prec);
            assert_abs_diff_eq!(out.2, expect.2, epsilon = prec);
        })
    }
}

use crate::proj::Proj;
use crate::transform::{transform, Transform};
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
