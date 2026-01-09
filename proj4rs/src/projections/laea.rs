//!
//! Lambert Azimuthal Equal Area
//!
//! ref: <https://proj.org/operations/projections/laea.html>
//!
//! laea:  "Lambert Azimuthal Equal Area" "\n\tAzi, Sph&Ell";
//!
use crate::errors::{Error, Result};
use crate::math::{
    authlat, authset,
    consts::{EPS_10, FRAC_PI_2, FRAC_PI_4},
    qsfn,
};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { laea }

#[derive(Debug, Clone)]
pub(crate) enum Projection {
    El(EProj),
    Sp(SProj),
}

impl Projection {
    pub fn laea(p: &mut ProjData, _: &ParamList) -> Result<Self> {
        Ok(if p.ellps.is_ellipsoid() {
            Self::El(EProj::new(p))
        } else {
            Self::Sp(SProj::new(p))
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::El(p) => p.forward(lam, phi, z),
            Self::Sp(p) => p.forward(lam, phi, z),
        }
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match self {
            Self::El(p) => p.inverse(x, y, z),
            Self::Sp(p) => p.inverse(x, y, z),
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }
}

//
// Ellipsoid
//

#[allow(non_camel_case_types)]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone)]
enum EMode {
    N_POLE,
    S_POLE,
    EQUIT {
        dd: f64,
        rq: f64,
        xmf: f64,
        ymf: f64,
    },
    OBLIQ {
        dd: f64,
        rq: f64,
        sinb1: f64,
        cosb1: f64,
        xmf: f64,
        ymf: f64,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct EProj {
    phi0: f64,
    e: f64,
    one_es: f64,
    qp: f64,
    apa: (f64, f64, f64),
    mode: EMode,
}

impl EProj {
    fn new(p: &mut ProjData) -> Self {
        let (e, one_es) = (p.ellps.e, p.ellps.one_es);
        let qp = qsfn(1., e, one_es);
        let apa = authset(p.ellps.es);
        let phi0 = p.phi0;

        use EMode::*;

        let t = phi0.abs();
        let mode = if (t - FRAC_PI_2).abs() < EPS_10 {
            if phi0 < 0. { S_POLE } else { N_POLE }
        } else if t.abs() < EPS_10 {
            let rq = (0.5 * qp).sqrt();
            EQUIT {
                rq,
                dd: 1. / rq,
                xmf: 1.,
                ymf: 0.5 * qp,
            }
        } else {
            let (sinphi, cosphi) = phi0.sin_cos();
            let rq = (0.5 * qp).sqrt();
            let sinb1 = qsfn(sinphi, e, one_es) / qp;
            let cosb1 = (1. - sinb1 * sinb1).sqrt();
            let dd = cosphi / ((1. - p.ellps.es * sinphi * sinphi).sqrt() * rq * cosb1);
            OBLIQ {
                rq,
                dd,
                sinb1,
                cosb1,
                xmf: rq * dd,
                ymf: rq / dd,
            }
        };

        Self {
            phi0,
            e,
            one_es,
            qp,
            apa,
            mode,
        }
    }

    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinlam, coslam) = lam.sin_cos();
        let q = qsfn(phi.sin(), self.e, self.one_es);

        use EMode::*;

        let (x, y) = match self.mode {
            OBLIQ {
                sinb1,
                cosb1,
                xmf,
                ymf,
                ..
            } => {
                let sinb = q / self.qp;
                let cosb = (1. - sinb * sinb).sqrt();
                let mut b = 1. + sinb1 * sinb + cosb1 * cosb * coslam;
                if b.abs() < EPS_10 {
                    return Err(Error::ToleranceConditionError);
                }
                b = (2. / b).sqrt();
                (
                    xmf * b * cosb * sinlam,
                    ymf * b * (cosb1 * sinb - sinb1 * cosb * coslam),
                )
            }
            EQUIT { xmf, ymf, .. } => {
                let sinb = q / self.qp;
                let cosb = (1. - sinb * sinb).sqrt();
                let mut b = 1. + cosb * coslam;
                if b.abs() < EPS_10 {
                    return Err(Error::ToleranceConditionError);
                }
                b = (2. / b).sqrt();
                (xmf * b * cosb * sinlam, ymf * b * sinb)
            }
            N_POLE => {
                if (FRAC_PI_2 + phi).abs() < EPS_10 {
                    return Err(Error::ToleranceConditionError);
                }
                let q = self.qp - q;
                if q >= 0. {
                    let b = q.sqrt();
                    (b * sinlam, -b * coslam)
                } else {
                    (0., 0.)
                }
            }
            S_POLE => {
                if (phi - FRAC_PI_2) < EPS_10 {
                    return Err(Error::ToleranceConditionError);
                }
                let q = self.qp + q;
                if q >= 0. {
                    let b = q.sqrt();
                    (b * sinlam, b * coslam)
                } else {
                    (0., 0.)
                }
            }
        };

        Ok((x, y, z))
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        use EMode::*;

        let (ab, xx, yy) = match self.mode {
            EQUIT { dd, rq, .. } => {
                let (x, y) = (x / dd, y * dd);
                let rho = x.hypot(y);
                if rho < EPS_10 {
                    return Ok((0., self.phi0, z));
                }
                let (sce, cce) = (2. * (0.5 * rho / rq).asin()).sin_cos();
                (y * sce / rho, x * sce, rho * cce)
            }
            OBLIQ {
                dd,
                rq,
                sinb1,
                cosb1,
                ..
            } => {
                let (x, y) = (x / dd, y * dd);
                let rho = x.hypot(y);
                if rho < EPS_10 {
                    return Ok((0., self.phi0, z));
                }
                let (sce, cce) = (2. * (0.5 * rho / rq).asin()).sin_cos();
                (
                    cce * sinb1 + y * sce * cosb1 / rho,
                    x * sce,
                    rho * cosb1 * cce - y * sinb1 * sce,
                )
            }
            N_POLE | S_POLE => {
                let q = x * x + y * y;
                if q == 0. {
                    return Ok((0., self.phi0, z));
                }
                let ab = 1. - q / self.qp;
                if matches!(self.mode, N_POLE) {
                    (ab, x, -y)
                } else {
                    (-ab, x, y)
                }
            }
        };
        Ok((xx.atan2(yy), authlat(ab.asin(), self.apa), z))
    }
}

//
// Sphere
//

#[allow(non_camel_case_types)]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone)]
enum SMode {
    N_POLE,
    S_POLE,
    EQUIT,
    OBLIQ { sinb1: f64, cosb1: f64 },
}

#[derive(Debug, Clone)]
pub(crate) struct SProj {
    phi0: f64,
    mode: SMode,
}

impl SProj {
    fn new(p: &mut ProjData) -> Self {
        let phi0 = p.phi0;
        let t = phi0.abs();

        use SMode::*;

        let mode = if (t - FRAC_PI_2).abs() < EPS_10 {
            if p.phi0 < 0. { S_POLE } else { N_POLE }
        } else if t.abs() < EPS_10 {
            EQUIT
        } else {
            let (sinb1, cosb1) = p.phi0.sin_cos();
            OBLIQ { sinb1, cosb1 }
        };
        Self { phi0, mode }
    }

    fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (sinphi, cosphi) = phi.sin_cos();
        let coslam = lam.cos();

        use SMode::*;

        match self.mode {
            EQUIT => {
                let mut y = 1. + cosphi * coslam;
                if y < EPS_10 {
                    Err(Error::ToleranceConditionError)
                } else {
                    y = (2. / y).sqrt();
                    Ok((y * cosphi * lam.sin(), y * sinphi, z))
                }
            }
            OBLIQ { sinb1, cosb1, .. } => {
                let mut y = 1. + sinb1 * sinphi + cosb1 * cosphi * coslam;
                if y < EPS_10 {
                    Err(Error::ToleranceConditionError)
                } else {
                    y = (2. / y).sqrt();
                    Ok((
                        y * cosphi * lam.sin(),
                        y * cosb1 * sinphi - sinb1 * cosphi * coslam,
                        z,
                    ))
                }
            }
            N_POLE | S_POLE => {
                if (phi + self.phi0).abs() < EPS_10 {
                    Err(Error::ToleranceConditionError)
                } else {
                    let mut y = FRAC_PI_4 - phi * 0.5;
                    if matches!(self.mode, N_POLE) {
                        y = 2. * y.sin();
                        Ok((y * lam.sin(), y * -coslam, z))
                    } else {
                        y = 2. * y.cos();
                        Ok((y * lam.sin(), y * coslam, z))
                    }
                }
            }
        }
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let rh = x.hypot(y);
        let mut phi = rh * 0.5;
        if phi > 1. {
            return Err(Error::ToleranceConditionError);
        }

        let lam;

        use SMode::*;

        phi = 2. * phi.asin();
        match self.mode {
            EQUIT => {
                let (sinz, cosz) = phi.sin_cos();
                phi = if rh <= EPS_10 {
                    0.
                } else {
                    (y * sinz / rh).asin()
                };
                let yy = cosz * rh;
                lam = if yy == 0. { 0. } else { (x * sinz).atan2(yy) };
            }
            OBLIQ { sinb1, cosb1 } => {
                let (sinz, cosz) = phi.sin_cos();
                phi = if rh <= EPS_10 {
                    self.phi0
                } else {
                    (cosz * sinb1 + y * sinz * cosb1 / rh).asin()
                };
                let yy = (cosz - phi.sin() * sinb1) * rh;
                lam = if yy == 0. {
                    0.
                } else {
                    (x * sinz * cosb1).atan2(yy)
                };
            }
            N_POLE => {
                phi = FRAC_PI_2 - phi;
                lam = x.atan2(-y);
            }
            S_POLE => {
                phi -= FRAC_PI_2;
                lam = x.atan2(y);
            }
        }
        Ok((lam, phi, z))
    }
}

#[cfg(test)]
mod tests {
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_laea_el() {
        let p = Proj::from_proj_string("+proj=laea +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222602.471450095181, 110589.82722441027, 0.)),
            (
                (2., -1., 0.),
                (222602.471450095181, -110589.827224408786, 0.),
            ),
            (
                (-2., 1., 0.),
                (-222602.471450095181, 110589.82722441027, 0.),
            ),
            (
                (-2., -1., 0.),
                (-222602.471450095181, -110589.827224408786, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, 1.0e-8);
    }

    #[test]
    fn proj_laea_sp() {
        let p = Proj::from_proj_string("+proj=laea +a=6400000").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (223365.281370124663, 111716.668072915665, 0.)),
            (
                (2., -1., 0.),
                (223365.281370124663, -111716.668072915665, 0.),
            ),
            (
                (-2., 1., 0.),
                (-223365.281370124663, 111716.668072915665, 0.),
            ),
            (
                (-2., -1., 0.),
                (-223365.281370124663, -111716.668072915665, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn test_epsg3035() {
        // cf https://github.com/3liz/proj4rs/issues/18

        let p = Proj::from_proj_string(
            "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80",
        )
        .unwrap();

        println!("{:#?}", p.projection());

        let inputs = [(
            (15.4213696, 47.0766716, 0.),
            (4732659.007426266, 2677630.7269610995, 0.),
        )];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, 1.0e-7);
    }
}
