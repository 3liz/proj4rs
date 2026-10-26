//!
//! Laborde
//!
//! ref: <https://proj.org/operations/projections/labrd.html>
//!
//! labrd: "Laborde" "\n\tCyl, Sph\n\tSpecial for Madagascar\n\tlat_0=";
//!
use crate::errors::{Error, Result};
use crate::math::consts::FRAC_PI_4;
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { labrd }

const EPS: f64 = 1.0e-10;

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    e: f64,
    one_es: f64,
    k0: f64,
    phi0: f64,
    k_rg: f64,
    p0s: f64,
    a: f64,
    c: f64,
    ca: f64,
    cb: f64,
    cc: f64,
    cd: f64,
}

impl Projection {
    pub fn labrd(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        if p.phi0 == 0. {
            return Err(Error::InvalidParameterValue(
                "Invalid value for lat_0: lat_0 should be different from 0",
            ));
        }

        let az = params.try_angular_value("azi")?.unwrap_or(0.);
        let el = &p.ellps;

        let sinp = p.phi0.sin();
        let t = 1. - el.es * sinp * sinp;
        let n = 1. / t.sqrt();
        let r = el.one_es * n / t;
        let k_rg = p.k0 * (n * r).sqrt();
        let p0s = ((r / n).sqrt() * p.phi0.tan()).atan();
        let a = sinp / p0s.sin();
        let t = el.e * sinp;
        let c = 0.5 * el.e * a * ((1. + t) / (1. - t)).ln()
            - a * (FRAC_PI_4 + 0.5 * p.phi0).tan().ln()
            + (FRAC_PI_4 + 0.5 * p0s).tan().ln();
        let t = az + az;
        let cb = 1. / (12. * k_rg * k_rg);
        let ca = (1. - t.cos()) * cb;
        let cb = cb * t.sin();

        Ok(Self {
            e: el.e,
            one_es: el.one_es,
            k0: p.k0,
            phi0: p.phi0,
            k_rg,
            p0s,
            a,
            c,
            ca,
            cb,
            cc: 3. * (ca * ca - cb * cb),
            cd: 6. * ca * cb,
        })
    }

    #[inline(always)]
    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let v1 = self.a * (FRAC_PI_4 + 0.5 * phi).tan().ln();
        let t = self.e * phi.sin();
        let v2 = 0.5 * self.e * self.a * ((1. + t) / (1. - t)).ln();
        let ps = 2. * ((v1 - v2 + self.c).exp().atan() - FRAC_PI_4);
        let i1 = ps - self.p0s;
        let (sinps, cosps) = ps.sin_cos();
        let cosps2 = cosps * cosps;
        let sinps2 = sinps * sinps;
        let i4 = self.a * cosps;
        let i2 = 0.5 * self.a * i4 * sinps;
        let i3 = i2 * self.a * self.a * (5. * cosps2 - sinps2) / 12.;
        let i6 = i4 * self.a * self.a;
        let i5 = i6 * (cosps2 - sinps2) / 6.;
        let i6 = i6
            * (self.a * self.a * (5. * cosps2 * cosps2 + sinps2 * (sinps2 - 18. * cosps2)) / 120.);
        let t = lam * lam;
        let x = self.k_rg * lam * (i4 + t * (i5 + t * i6));
        let y = self.k_rg * (i1 + t * (i2 + t * i3));
        let x2 = x * x;
        let y2 = y * y;
        let v1 = 3. * x * y2 - x * x2;
        let v2 = y * y2 - 3. * x2 * y;
        Ok((
            x + self.ca * v1 + self.cb * v2,
            y + self.ca * v2 - self.cb * v1,
            z,
        ))
    }

    #[inline(always)]
    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let x2 = x * x;
        let y2 = y * y;
        let v1 = 3. * x * y2 - x * x2;
        let v2 = y * y2 - 3. * x2 * y;
        let v3 = x * (5. * y2 * y2 + x2 * (-10. * y2 + x2));
        let v4 = y * (5. * x2 * x2 + y2 * (-10. * x2 + y2));
        let x = x - self.ca * v1 - self.cb * v2 + self.cc * v3 + self.cd * v4;
        let y = y + self.cb * v1 - self.ca * v2 - self.cd * v3 + self.cc * v4;
        let ps = self.p0s + y / self.k_rg;
        let mut pe = ps + self.phi0 - self.p0s;

        for _ in 0..20 {
            let v1 = self.a * (FRAC_PI_4 + 0.5 * pe).tan().ln();
            let tpe = self.e * pe.sin();
            let v2 = 0.5 * self.e * self.a * ((1. + tpe) / (1. - tpe)).ln();
            let t = ps - 2. * ((v1 - v2 + self.c).exp().atan() - FRAC_PI_4);
            pe += t;
            if t.abs() < EPS {
                break;
            }
        }

        let t = self.e * pe.sin();
        let t = 1. - t * t;
        let re = self.one_es / (t * t.sqrt());
        let t = ps.tan();
        let t2 = t * t;
        let s = self.k_rg * self.k_rg;
        let d = re * self.k0 * self.k_rg;
        let i7 = t / (2. * d);
        let i8 = t * (5. + 3. * t2) / (24. * d * s);
        let d = ps.cos() * self.k_rg * self.a;
        let i9 = 1. / d;
        let d = d * s;
        let i10 = (1. + 2. * t2) / (6. * d);
        let i11 = (5. + t2 * (28. + 24. * t2)) / (120. * d * s);
        let x2 = x * x;
        Ok((
            x * (i9 + x2 * (-i10 + x2 * i11)),
            pe + x2 * (-i7 + i8 * x2),
            z,
        ))
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

    // Test values from proj builtins.gie
    #[test]
    fn proj_labrd() {
        let p = Proj::from_proj_string("+proj=labrd +ellps=GRS80 +lon_0=0.5 +lat_0=2").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (166973.166090228, -110536.912730266, 0.)),
            ((2., -1., 0.), (166973.168287157, -331761.993650884, 0.)),
            ((-2., 1., 0.), (-278345.500519976, -110469.032642032, 0.)),
            ((-2., -1., 0.), (-278345.504185270, -331829.870790275, 0.)),
        ];
        test_proj_forward(&p, &inputs, 1e-4);

        let inputs = [
            ((0.501797719, 2.000904357, 0.), (200., 100., 0.)),
            ((0.501797717, 1.999095641, 0.), (200., -100., 0.)),
            ((0.498202281, 2.000904357, 0.), (-200., 100., 0.)),
            ((0.498202283, 1.999095641, 0.), (-200., -100., 0.)),
        ];
        test_proj_inverse(&p, &inputs, 1e-9);
    }

    #[test]
    fn proj_labrd_lat_0() {
        assert!(Proj::from_proj_string("+proj=labrd +ellps=GRS80 +lat_0=0").is_err());
    }
}
