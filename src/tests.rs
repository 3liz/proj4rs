//!
//! Unit tests
//!

pub(crate) mod utils {
    use crate::errors::Result;
    use crate::proj::{Proj, ProjData};
    use approx::assert_abs_diff_eq;

    pub fn scale(d: &ProjData, xyz: (f64, f64, f64)) -> (f64, f64, f64) {
        (xyz.0 * d.ellps.a + d.x0 , xyz.1 * d.ellps.a + d.y0, xyz.2)
    }

    pub fn descale(d: &ProjData, xyz: (f64, f64, f64)) -> (f64, f64, f64) {
        ((xyz.0 - d.x0) * d.ellps.ra, (xyz.1 - d.y0) * d.ellps.ra, xyz.2)
    }

    pub fn to_deg(lam: f64, phi: f64, z: f64) -> (f64, f64, f64) {
        (lam.to_degrees(), phi.to_degrees(), z)
    }

    pub fn to_rad(lpz: (f64, f64, f64)) -> (f64, f64, f64) {
        (lpz.0.to_radians(), lpz.1.to_radians(), lpz.2)
    }

    pub fn test_proj_forward(p: &Proj, inputs: &[((f64, f64, f64), (f64, f64, f64))], prec: f64) {
        let d = p.data();
        inputs.iter().for_each(|(input, expect)| {
            let (lam, phi, z) = to_rad(*input);
            let out = scale(
                d,
                p.projection().forward(lam - d.lam0, phi, z).unwrap(),
            );
            assert_abs_diff_eq!(out.0, expect.0, epsilon = prec);
            assert_abs_diff_eq!(out.1, expect.1, epsilon = prec);
            assert_abs_diff_eq!(out.2, expect.2, epsilon = prec);
        })
    }

    pub fn test_proj_inverse(p: &Proj, inputs: &[((f64, f64, f64), (f64, f64, f64))], prec: f64) {
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
