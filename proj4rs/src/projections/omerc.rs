//!
//! Oblique Mercator
//!
//! ref: <https://proj.org/en/stable/operations/projections/omerc.html>
//! 
//! ```text
//! omerc: "Oblique Mercator" "\n\tCyl, Sph&Ell no_rot\n\t"
//!     "alpha= [gamma=] [no_off] lonc= or\n\t lon_1= lat_1= lon_2= lat_2=";
//! ```
//!
//! Ported from PROJ `src/projections/omerc.cpp`.
//!
//! Copyright (c) 2003, 2006   Gerald I. Evenden
//!
//! Permission is hereby granted, free of charge, to any person obtaining
//! a copy of this software and associated documentation files (the
//! "Software"), to deal in the Software without restriction, including
//! without limitation the rights to use, copy, modify, merge, publish,
//! distribute, sublicense, and/or sell copies of the Software, and to
//! permit persons to whom the Software is furnished to do so, subject to
//! the following conditions:
//!
//! The above copyright notice and this permission notice shall be
//! included in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//!
use crate::errors::{Error, Result};
use crate::math::{
    self, aasin, adjlon,
    consts::{EPS_10, FRAC_PI_2, FRAC_PI_4, PI, TAU},
    tsfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { omerc }

const TOL: f64 = 1.0e-7;

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub(crate) struct Projection {
    e: f64,
    A: f64,
    B: f64,
    E: f64,
    ArB: f64,
    BrA: f64,
    rB: f64,
    singam: f64,
    cosgam: f64,
    sinrot: f64,
    cosrot: f64,
    v_pole_n: f64,
    v_pole_s: f64,
    u_0: f64,
    no_rot: bool,
}

#[allow(non_snake_case)]
impl Projection {
    pub fn omerc(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        let (e, es, one_es) = (p.ellps.e, p.ellps.es, p.ellps.one_es);
        let phi0 = p.phi0;

        let no_rot = params.check_option("no_rot")?;
        let alpha_p = params.try_angular_value("alpha")?;
        let gamma_p = params.try_angular_value("gamma")?;

        let (mut lamc, mut no_off) = (0., false);
        let (mut lam1, mut phi1, mut lam2, mut phi2) = (0., 0., 0., 0.);
        if alpha_p.is_some() || gamma_p.is_some() {
            lamc = params.try_angular_value("lonc")?.unwrap_or(0.);
            no_off = params.check_option("no_off")? || params.check_option("no_uoff")?;
        } else {
            lam1 = params.try_angular_value("lon_1")?.unwrap_or(0.);
            phi1 = params.try_angular_value("lat_1")?.unwrap_or(0.);
            lam2 = params.try_angular_value("lon_2")?.unwrap_or(0.);
            phi2 = params.try_angular_value("lat_2")?.unwrap_or(0.);

            if phi1.abs() > FRAC_PI_2 - TOL {
                return Err(Error::InvalidParameterValue("|lat_1| should be < 90°"));
            }
            if phi2.abs() > FRAC_PI_2 - TOL {
                return Err(Error::InvalidParameterValue("|lat_2| should be < 90°"));
            }
            if (phi1 - phi2).abs() <= TOL {
                return Err(Error::InvalidParameterValue(
                    "lat_1 should be different from lat_2",
                ));
            }
            if phi1.abs() <= TOL {
                return Err(Error::InvalidParameterValue(
                    "lat_1 should be different from 0",
                ));
            }
        }
        if (phi0.abs() - FRAC_PI_2).abs() <= TOL {
            return Err(Error::InvalidParameterValue("|lat_0| should be < 90°"));
        }

        let com = one_es.sqrt();
        let (A, B, E, D, mut F);
        if phi0.abs() > EPS_10 {
            let (sinph0, cosph0) = phi0.sin_cos();
            let con = 1. - es * sinph0 * sinph0;
            let b = cosph0 * cosph0;
            B = (1. + es * b * b / one_es).sqrt();
            A = B * p.k0 * com / con;
            D = B * com / (cosph0 * con.sqrt());
            F = D * D - 1.;
            if F <= 0. {
                F = 0.;
            } else {
                F = F.sqrt();
                if phi0 < 0. {
                    F = -F;
                }
            }
            F += D;
            E = F * tsfn(phi0, sinph0, e).powf(B);
        } else {
            B = 1. / com;
            A = p.k0;
            D = 1.;
            F = 1.;
            E = 1.;
        }

        let (gamma0, gamma, alpha_c);
        if alpha_p.is_some() || gamma_p.is_some() {
            if let Some(alpha) = alpha_p {
                gamma0 = aasin(alpha.sin() / D)?;
                gamma = gamma_p.unwrap_or(alpha);
                alpha_c = alpha;
            } else {
                let g = gamma_p.unwrap();
                gamma0 = g;
                gamma = g;
                // For a sphere, |gamma| must be <= 90 - |lat_0|
                alpha_c = aasin(D * gamma0.sin())
                    .map_err(|_| Error::InvalidParameterValue("|gamma| too large for lat_0"))?;
            }
            p.lam0 = lamc - aasin(0.5 * (F - 1. / F) * gamma0.tan())? / B;
        } else {
            let H = tsfn(phi1, phi1.sin(), e).powf(B);
            let L = tsfn(phi2, phi2.sin(), e).powf(B);
            F = E / H;
            let pp = (L - H) / (L + H);
            if pp == 0. {
                return Err(Error::InvalidParameterValue(
                    "Invalid value for eccentricity",
                ));
            }
            let J = E * E;
            let J = (J - L * H) / (J + L * H);
            let con = lam1 - lam2;
            if con < -PI {
                lam2 -= TAU;
            } else if con > PI {
                lam2 += TAU;
            }
            p.lam0 =
                adjlon(0.5 * (lam1 + lam2) - (J * (0.5 * B * (lam1 - lam2)).tan() / pp).atan() / B);
            let denom = F - 1. / F;
            if denom == 0. {
                return Err(Error::InvalidParameterValue(
                    "Invalid value for eccentricity",
                ));
            }
            gamma0 = (2. * (B * adjlon(lam1 - p.lam0)).sin() / denom).atan();
            alpha_c = aasin(D * gamma0.sin())?;
            gamma = alpha_c;
        }

        let rB = 1. / B;
        let ArB = A * rB;
        let u_0 = if no_off {
            0.
        } else {
            let u_0 = (ArB * ((D * D - 1.).sqrt() / alpha_c.cos()).atan()).abs();
            if phi0 < 0. { -u_0 } else { u_0 }
        };
        let (singam, cosgam) = gamma0.sin_cos();
        let (sinrot, cosrot) = gamma.sin_cos();
        let F = 0.5 * gamma0;

        Ok(Self {
            e,
            A,
            B,
            E,
            ArB,
            BrA: 1. / ArB,
            rB,
            singam,
            cosgam,
            sinrot,
            cosrot,
            v_pole_n: ArB * (FRAC_PI_4 - F).tan().ln(),
            v_pole_s: ArB * (FRAC_PI_4 + F).tan().ln(),
            u_0,
            no_rot,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (mut u, v);
        if (phi.abs() - FRAC_PI_2).abs() > EPS_10 {
            let W = self.E / tsfn(phi, phi.sin(), self.e).powf(self.B);
            let one_div_W = 1. / W;
            let S = 0.5 * (W - one_div_W);
            let T = 0.5 * (W + one_div_W);
            let V = (self.B * lam).sin();
            let U = (S * self.singam - V * self.cosgam) / T;
            if (U.abs() - 1.).abs() < EPS_10 {
                return Err(Error::CoordinateOutOfRange);
            }
            v = 0.5 * self.ArB * ((1. - U) / (1. + U)).ln();
            let temp = (self.B * lam).cos();
            u = if temp.abs() < TOL {
                self.A * lam
            } else {
                self.ArB * (S * self.cosgam + V * self.singam).atan2(temp)
            };
        } else {
            v = if phi > 0. {
                self.v_pole_n
            } else {
                self.v_pole_s
            };
            u = self.ArB * phi;
        }
        if self.no_rot {
            Ok((u, v, z))
        } else {
            u -= self.u_0;
            Ok((
                v * self.cosrot + u * self.sinrot,
                u * self.cosrot - v * self.sinrot,
                z,
            ))
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (u, v) = if self.no_rot {
            (x, y)
        } else {
            (
                y * self.cosrot + x * self.sinrot + self.u_0,
                x * self.cosrot - y * self.sinrot,
            )
        };
        let Qp = (-self.BrA * v).exp();
        if Qp == 0. {
            return Err(Error::CoordinateOutOfRange);
        }
        let Sp = 0.5 * (Qp - 1. / Qp);
        let Tp = 0.5 * (Qp + 1. / Qp);
        let Vp = (self.BrA * u).sin();
        let Up = (Vp * self.cosgam + Sp * self.singam) / Tp;
        if (Up.abs() - 1.).abs() < EPS_10 {
            Ok((0., FRAC_PI_2.copysign(Up), z))
        } else {
            let ts = (self.E / ((1. + Up) / (1. - Up)).sqrt()).powf(1. / self.B);
            let phi = math::phi2(ts, self.e)?;
            let lam = -self.rB * (Sp * self.cosgam - Vp * self.singam).atan2((self.BrA * u).cos());
            Ok((lam, phi, z))
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

//============
// Tests
//============

#[cfg(test)]
mod tests {
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    // Expected values computed with PROJ 9.8.1 (through pyproj 3.8.0).

    #[test]
    fn proj_omerc_alpha_gamma_no_uoff() {
        // EPSG:3375 (Malaysia RSO, variant A)
        let p = Proj::from_proj_string(concat!(
            "+proj=omerc +lat_0=4 +lonc=102.25 +alpha=323.0257964666666 +k=0.99984",
            " +x_0=804671 +y_0=0 +no_uoff +gamma=323.1301023611111 +ellps=GRS80"
        ))
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((102.25, 4., 0.), (472830.426071632, 442454.098683051, 0.)),
            ((101.5, 3.1, 0.), (389287.174191490, 343128.619662765, 0.)),
            ((104., 1.5, 0.), (667104.332086940, 165868.754500688, 0.)),
            ((100.5, 6.5, 0.), (279734.832387836, 719467.367720329, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        // EPSG:26931 (Alaska zone 1)
        let p = Proj::from_proj_string(concat!(
            "+proj=omerc +lat_0=57 +lonc=-133.6666666666667 +alpha=323.1301023611111 +k=0.9999",
            " +x_0=5000000 +y_0=-5000000 +no_uoff +gamma=323.1301023611111 +ellps=GRS80"
        ))
        .unwrap();

        let inputs = [
            (
                (-133.6666666666667, 57., 0.),
                (818676.734401053, 575097.688875169, 0.),
            ),
            ((-135., 58., 0.), (739841.600324485, 687227.401595449, 0.)),
            ((-131., 55.5, 0.), (987154.626888720, 411358.336246798, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_alpha_gamma() {
        // EPSG:29873 (Borneo RSO, variant B)
        let p = Proj::from_proj_string(concat!(
            "+proj=omerc +lat_0=4 +lonc=115 +alpha=53.31582047222222 +k=0.99984",
            " +x_0=590476.87 +y_0=442857.65 +gamma=53.13010236111111 +ellps=evrstSS"
        ))
        .unwrap();

        let inputs = [
            ((115., 4., 0.), (590476.87, 442857.65, 0.)),
            ((117., 6., 0.), (811179.866049288, 664991.329329078, 0.)),
            ((110., 1.5, 0.), (35001.128040784, 165849.219077694, 0.)),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_alpha() {
        let p =
            Proj::from_proj_string("+proj=omerc +lat_0=0 +lonc=0 +alpha=45 +ellps=WGS84").unwrap();
        let inputs = [
            ((0., 0., 0.), (0., 0., 0.)),
            ((10., -5., 0.), (1117856.956841016, -556763.900779827, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);

        let p = Proj::from_proj_string(
            "+proj=omerc +lat_0=-30 +lonc=140 +alpha=-60 +k=0.9996 +ellps=GRS80",
        )
        .unwrap();
        let inputs = [
            ((140., -30., 0.), (0., 0., 0.)),
            ((150., -35., 0.), (912497.393815875, -595944.407570422, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_gamma() {
        let p = Proj::from_proj_string("+proj=omerc +lat_0=40 +lonc=10 +gamma=20 +ellps=WGS84")
            .unwrap();
        let inputs = [
            ((10., 40., 0.), (0., 0., 0.)),
            ((12., 42., 0.), (139309.854872632, 241309.802316528, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_no_rot() {
        let p =
            Proj::from_proj_string("+proj=omerc +lat_0=40 +lonc=10 +alpha=20 +no_rot +ellps=WGS84")
                .unwrap();
        let inputs = [
            ((10., 40., 0.), (4640083.832007261, 0., 0.)),
            ((12., 42., 0.), (4907261.884510875, 79095.125455320, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_two_points() {
        let p = Proj::from_proj_string(
            "+proj=omerc +lat_1=45 +lon_1=-100 +lat_2=50 +lon_2=-80 +lat_0=47 +ellps=WGS84",
        )
        .unwrap();
        let inputs = [
            ((-90., 47., 0.), (298828.510535889, 7536.787031617, 0.)),
            ((-100., 45., 0.), (-478192.174780696, -203962.733704721, 0.)),
            ((-80., 50., 0.), (993675.885257817, 423831.381319845, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_sphere() {
        let p = Proj::from_proj_string("+proj=omerc +R=6370997 +lat_0=-42 +lonc=147 +alpha=30")
            .unwrap();
        let inputs = [
            ((147., -42., 0.), (0., 0., 0.)),
            ((150., -30., 0.), (294093.952081180, 1329914.798619571, 0.)),
            ((140., -50., 0.), (-501386.425950173, -911093.204344402, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }

    #[test]
    fn proj_omerc_invalid() {
        for def in [
            "+proj=omerc +lat_1=45 +lon_1=-100 +lat_2=45 +lon_2=-80 +lat_0=47 +ellps=WGS84",
            "+proj=omerc +lat_0=90 +lonc=0 +alpha=10 +ellps=WGS84",
            "+proj=omerc +lat_0=60 +lonc=0 +gamma=50 +ellps=WGS84",
        ] {
            assert!(Proj::from_proj_string(def).is_err(), "{def}");
        }
    }
}
