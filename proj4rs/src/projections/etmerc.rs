//!
//! Implement the Knud Poder/Karsten Engsager algorithm.
//!
//! This algorithm is called "exact" in comparison to the Evenden/John Snyder algorithm (tmerc)
//! (slightly) faster but less accurate. Note that etmerc does not provide spherical projection.
//!
//! This algorithm is used as a base for UTM projections.
//!
//! Reference: <https://proj.org/operations/projections/tmerc.html>
//!
//! etmerc: "Extended Transverse Mercator" "\n\tCyl, Sph\n\tlat_ts=(0)\nlat_0=(0)"
//! utm: "Universal Transverse Mercator (UTM)" "\n\tCyl, Sph\n\tzone= south"
//!
#![allow(non_snake_case)]

// Projection stub
super::projection! { etmerc, utm }

use crate::errors::{Error, Result};
use crate::math::{adjlon, asinh, consts::PI};
use crate::parameters::ParamList;
use crate::proj::ProjData;

const ETMERC_ORDER: usize = 6;

type Coeffs = [f64; ETMERC_ORDER];

#[inline]
fn gatg(c: &Coeffs, B: f64) -> f64 {
    let mut h2 = 0.;
    let cos_2B = 2. * (2. * B).cos();
    let h = c
        .iter()
        .copied()
        .reduce(|h1, p| {
            let h = -h2 + cos_2B * h1 + p;
            h2 = h1;
            h
        })
        .unwrap();
    B + h * (2. * B).sin()
}

// Complex Clenshaw summation
#[inline]
fn clens_cplx(a: &Coeffs, arg_r: f64, arg_i: f64) -> (f64, f64) {
    let (sin_arg_r, cos_arg_r) = arg_r.sin_cos();
    let sinh_arg_i = arg_i.sinh();
    let cosh_arg_i = arg_i.cosh();

    let mut r = 2. * cos_arg_r * cosh_arg_i;
    let mut i = -2. * sin_arg_r * sinh_arg_i;

    let (mut hi1, mut hi2) = (0., 0.);
    let (mut hr1, mut hr2) = (0., 0.);

    let mut hi = 0.;
    let hr = a
        .iter()
        .copied()
        .reduce(|hr, p| {
            hr2 = hr1;
            hi2 = hi1;
            hr1 = hr;
            hi1 = hi;
            hi = -hi2 + i * hr1 + r * hi1;
            -hr2 + r * hr1 - i * hi1 + p
        })
        .unwrap();

    r = sin_arg_r * cosh_arg_i;
    i = cos_arg_r * sinh_arg_i;
    (
        r * hr - i * hi, // R
        r * hi + i * hr, // I
    )
}

// Real Clenshaw summation
#[inline]
fn clens(a: &Coeffs, arg_r: f64) -> f64 {
    let cos_arg_r = arg_r.cos();
    let r = 2. * cos_arg_r;

    let (mut hr1, mut hr2) = (0., 0.);
    let hr = unsafe {
        a.iter()
            .copied()
            .reduce(|hr, p| {
                hr2 = hr1;
                hr1 = hr;
                -hr2 + r * hr1 + p
            })
            .unwrap_unchecked()
    };
    arg_r.sin() * hr
}

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    Qn: f64,     // Merid. quad., scaled to the projection
    Zb: f64,     // Radius vector in polar coord. systems
    cgb: Coeffs, // Constants for Gauss -> Geo lat
    cbg: Coeffs, // Constants for Geo lat -> Gauss
    utg: Coeffs, // Constants for transv. merc. -> geo
    gtu: Coeffs, // Constants for geo -> transv. merc.
}

#[rustfmt::skip]
impl Projection {
    pub fn etmerc(p: &mut ProjData, _params: &ParamList) -> Result<Self> {

        // We have flattening computed, use it !
        let f = p.ellps.f;

        if f == 0. {
            return Err(Error::EllipsoidRequired)
        }

        // third flattening
        let n = f / (2. - f);

        // COEF. OF TRIG SERIES GEO <-> GAUSS */
        // cgb := Gaussian -> Geodetic, KW p190 - 191 (61) - (62) */
        // cbg := Geodetic -> Gaussian, KW p186 - 187 (51) - (52) */
        // ETMERC_ORDER = 6th degree : Engsager and Poder: ICC2007 */
        let n2 = n*n;

        let mut cgb: Coeffs = [
            n * ( 2. + n*(-2./3.0 + n*(-2. + n*(116./45.0 + n*(26./45.0 + n*(-2854./675.0 )))))),
            n2 * (7./3.0 + n*(-8./5.0 + n*(-227./45.0 + n*(2704./315.0 + n*(2323./945.0))))),
            n2 * n * (56./15.0 + n*(-136./35.0 + n*(-1262./105.0 + n*(73814./2835.0)))),
            n2 * n2 * (4279./630.0 + n*(-332./35.0 + n*(-399572./14175.0))),
            n2 * n2 * n * (4174./315.0 + n*(-144838./6237.0)),
            n2 * n2 * n2 * (601676./22275.0),
        ];

        let mut cbg: Coeffs = [ 
            n * (-2. + n*( 2./3.0 + n*( 4./3.0  + n*(-82./45.0 + n*(32./45.0 + n*(4642./4725.0)))))),
            n2 * (5./3.0 + n*(-16./15.0 + n*(-13./9.0  + n*(904./315.0 + n*(-1522./945.0))))),
            n2 * n * (-26./15.0 + n*(34./21.0 + n*(8./5.0 + n*(-12686./2835.0)))),
            n2 * n2 * (1237./630.0 + n*(-12./5.0 + n*(-24832./14175.0))),
            n2 * n2 * n * (-734./315.0 + n*(109598./31185.0)),
            n2 * n2 * n2 * (444337./155925.0),
        ];

        // Coefficients are used backward so reverse them now 
        cgb.reverse();
        cbg.reverse();

        // Constants of the projections 
        // Transverse Mercator (UTM, ITM, etc)
        
        // Norm. mer. quad, K&W p.50 (96), p.19 (38b), p.5 (2)
        let Qn = p.k0/(1. + n) * (1. + n2*(1./4.0 + n2*(1./64.0 + n2/256.0)));
    
        // coef of trig series 
        // utg := ell. N, E -> sph. N, E,  KW p194 (65)
        // gtu := sph. N, E -> ell. N, E,  KW p196 (69)
        let mut utg: Coeffs = [
            n * (-0.5  + n*(2./3.0 + n*(-37./96.0 + n*(1./360.0 + n*(81./512.0 + n*(-96199./604800.0)))))),
            n2 * (-1./48.0 + n*(-1./15.0 + n*(437./1440.0 + n*(-46./105.0 + n*(1118711./3870720.0))))),
            n2 * n *(-17./480.0 + n*( 37./840.0 + n*(209./4480.0 + n*(-5569./90720.0)))),
            n2 * n2 *(-4397./161280.0 + n*(11./504.0 + n*( 830251./7257600.0))),
            n2 * n2 * n * (-4583./161280.0 + n*(108847./3991680.0)),
            n2 * n2 *n2 *(-20648693./638668800.0),
        ];

        let mut gtu: Coeffs = [
            n *( 0.5 + n*(-2./3.0 + n*(5./16.0 + n*(41./180.0 + n*(-127./288.0 + n*(7891./37800.0)))))),
            n2 * (13./48.0 + n*(-3./5.0 + n*(557./1440.0 + n*(281./630.0 + n*(-1983433./1935360.0))))),
            n2 * n * (61./240.0 + n*(-103./140.0 + n*(15061./26880.0 + n*(167603./181440.0)))),
            n2 * n2 * (49561./161280.0 + n*(-179./168.0 + n*(6601661./7257600.0))),
            n2 * n2 * n * (34729./80640.0 + n*(-3418889./1995840.0)),
            n2 * n2 * n2 * (212378941./319334400.0),
        ];

        // Coefficients are used backward so reverse them now 
        utg.reverse();
        gtu.reverse();

        // Gaussian latitude value of the origin latitude
        let z = gatg(&cbg, p.phi0);

        // Origin northing minus true northing at the origin latitude 
        // i.e. true northing = N - Zb
        let Zb  = - Qn *(z + clens(&gtu, 2.*z));

        Ok(Self {
            Qn,
            Zb,
            cgb,
            cbg,
            utg,
            gtu,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {

        let (mut Cn, mut Ce) = (phi, lam);

        // ell. LAT, LNG -> Gaussian LAT, LNG
        Cn = gatg(&self.cbg, Cn);

        // Gaussian LAT, LNG -> compl. sph. LAT
        let (sin_Cn, cos_Cn) = Cn.sin_cos();
        let (sin_Ce, cos_Ce) = Ce.sin_cos();
        
        Cn = sin_Cn.atan2(cos_Ce*cos_Cn);
        Ce = (sin_Ce*cos_Cn).atan2(sin_Cn.hypot(cos_Cn*cos_Ce));

        // compl. sph. N, E -> ell. norm. N, E
        Ce  = asinh(Ce.tan());
        let (dCn, dCe) = clens_cplx(&self.gtu, 2.*Cn, 2.*Ce);
        Cn += dCn;
        Ce += dCe;

        if Ce.abs() <= 2.623395162778 {
            Ok((
                self.Qn * Ce,            // Easting
                self.Qn * Cn + self.Zb,  // Northing 
                z,
            ))
        } else {
            Err(Error::CoordTransOutsideProjectionDomain)
        }
    }

    pub fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (mut Cn, mut Ce) = (y, x);
        
        // normalize N, E
        Cn = (Cn - self.Zb)/self.Qn;
        Ce /= self.Qn;

        if Ce.abs() <= 2.623395162778 { // 150 degrees
            // norm. N, E -> compl. sph. LAT, LNG
            let (dCn, dCe) = clens_cplx(&self.utg, 2.*Cn, 2.*Ce);
            Cn += dCn;
            Ce += dCe;
            Ce = Ce.sinh().atan(); 
            // compl. sph. LAT -> Gaussian LAT, LNG
            let (sin_Cn, cos_Cn) = Cn.sin_cos();
            let (sin_Ce, cos_Ce) = Ce.sin_cos();

            Ce = sin_Ce.atan2(cos_Ce*cos_Cn);
            Cn = (sin_Cn*cos_Ce).atan2(sin_Ce.hypot(cos_Ce*cos_Cn));
            // Gaussian LAT, LNG -> ell. LAT, LNG
            Ok((
                Ce,
                gatg(&self.cgb, Cn),
                z,
            ))
        } else {
            Err(Error::InverseProjectionFailure)
        }
    }

    pub const fn has_inverse() -> bool {
        true
    }

    pub const fn has_forward() -> bool {
        true
    }

    //-------------------
    // UTM
    //------------------
    pub fn utm(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        p.x0 = 500_000.;
        p.y0 = if params.check_option("south")? {
            10_000_000.
        } else {
            0.
        };

        let zone = params.try_value::<i32>("zone").and_then(|zone| match zone {
            Some(zone) => {
                if (1..=60).contains(&zone) {
                    Ok(zone as f64)
                } else {
                    Err(Error::InvalidUtmZone)
                }
            }
            None => {
                // nearest central meridian input
                let zone = ((adjlon(p.lam0) + PI) * 30. / PI).floor().round();
                if (1. ..=60.).contains(&zone) {
                    Ok(zone)
                } else {
                    Err(Error::InvalidUtmZone)
                }
            }
        })?;

        p.lam0 = ((zone - 1.) + 0.5) * PI / 30. - PI;
        p.k0 = 0.9996;
        p.phi0 = 0.;

        Self::etmerc(p, params)
    }

}

#[cfg(test)]
mod tests {
    use crate::math::consts::EPS_10;
    use crate::proj::Proj;
    use crate::tests::utils::{test_proj_forward, test_proj_inverse};

    #[test]
    fn proj_etmerc_etmerc() {
        let p = Proj::from_proj_string("+proj=etmerc +ellps=GRS80").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            ((2., 1., 0.), (222650.79679758527, 110642.22941193319, 0.)),
            ((2., -1., 0.), (222650.79679758527, -110642.22941193319, 0.)),
            ((-2., 1., 0.), (-222650.79679758527, 110642.22941193319, 0.)),
            (
                (-2., -1., 0.),
                (-222650.79679758527, -110642.22941193319, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }

    #[test]
    fn proj_etmerc_utm() {
        let p = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();

        println!("{:#?}", p.projection());
        println!("{:#?}", p.data());

        let inputs = [
            ((2., 1., 0.), (1057002.4054912976, 110955.14117594929, 0.)),
            ((2., -1., 0.), (1057002.4054912976, -110955.1411759492, 0.)),
            ((-2., 1., 0.), (611263.8122789060, 110547.10569680421, 0.)),
            ((-2., -1., 0.), (611263.8122789060, -110547.10569680421, 0.)),
        ];

        test_proj_forward(&p, &inputs, EPS_10);
        test_proj_inverse(&p, &inputs, EPS_10);
    }
}
