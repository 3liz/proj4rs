//!
//! Project:  PROJ
//! Purpose:  Implementation of the krovak (Krovak) projection.
//!           Definition: http://www.ihsenergy.com/epsg/guid7.html#1.4.3
//! Author:   Thomas Flemming, tf@ttqv.com
//!
//!
//! Copyright (c) 2001, Thomas Flemming, tf@ttqv.com
//!
//! Permission is hereby granted, free of charge, to any person obtaining a
//! copy of this software and associated documentation files (the "Software"),
//! to deal in the Software without restriction, including without limitation
//! the rights to use, copy, modify, merge, publish, distribute, sublicense,
//! and/or sell copies of the Software, and to permit persons to whom the
//! Software is furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included
//! in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//! BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//! ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//! SOFTWARE.
//!
//!
//! A description of the (forward) projection is found in:
//!
//!```text
//!      Bohuslav Veverka,
//!
//!      KROVAK’S PROJECTION AND ITS USE FOR THE
//!      CZECH REPUBLIC AND THE SLOVAK REPUBLIC,
//!
//!      50 years of the Research Institute of
//!      and the Slovak Republic Geodesy, Topography and Cartography
//!```
//!
//! which can be found via the Wayback Machine:
//!
//!```text
//!      https://web.archive.org/web/20150216143806/https://www.vugtk.cz/odis/sborniky/sb2005/Sbornik_50_let_VUGTK/Part_1-Scientific_Contribution/16-Veverka.pdf
//!```
//!
//! Further info, including the inverse projection, is given by EPSG:
//!
//!```text
//!      Guidance Note 7 part 2
//!      Coordinate Conversions and Transformations including Formulas
//!
//!      http://www.iogp.org/pubs/373-07-2.pdf
//!```
//!
//! Variable names in this file mostly follows what is used in the
//! paper by Veverka.
//!
//! According to EPSG the full Krovak projection method should have
//! the following parameters.  Within PROJ the azimuth, and pseudo
//! standard parallel are hardcoded in the algorithm and can't be
//! altered from outside. The others all have defaults to match the
//! common usage with Krovak projection.
//!
//!```text
//!      lat_0 = latitude of centre of the projection
//!
//!      lon_0 = longitude of centre of the projection
//!
//!      ** = azimuth (true) of the centre line passing through the
//!           centre of the projection
//!
//!      ** = latitude of pseudo standard parallel
//!
//!      k  = scale factor on the pseudo standard parallel
//!
//!      x_0 = False Easting of the centre of the projection at the
//!            apex of the cone
//!
//!      y_0 = False Northing of the centre of the projection at
//!            the apex of the cone
//!```

use crate::ellps::{Ellipsoid, Shape};
use crate::errors::{Error, Result};
use crate::math::consts::{FRAC_PI_2, FRAC_PI_4};
use crate::parameters::ParamList;
use crate::proj::ProjData;

// Projection stub
super::projection! { krovak }

const EPS: f64 = 1.0e-15;
const UQ: f64 = 1.04216856380474; // DU(2, 59, 42, 42.69689)
const S0: f64 = 1.37008346281555; // Latitude of pseudo standard parallel 78deg 30'00" N

const MAX_ITER: usize = 100;

#[derive(Debug, Clone)]
pub(crate) struct Projection {
    e: f64,
    xyfact: (f64, f64),
    alpha: f64,
    k: f64,
    n: f64,
    rho0: f64,
    ad: f64,
    easting_northing: bool,
}

impl Projection {
    pub fn krovak(p: &mut ProjData, params: &ParamList) -> Result<Self> {
        // Bessel as fixed ellipsoid

        // NOTE: if we use the BESSEL definition from inverse
        // flattening we have a small difference (about 1.e-7 precision)
        // from output from Proj
        //p.ellps = Ellipsoid::try_from_ellipsoid(&BESSEL)?;
        p.ellps = Ellipsoid::calc_ellipsoid_params(6377397.155, Shape::SP_es(0.006674372230614))?;

        // If latitude of projection center is not set, use 49d30'N
        if params.get("lat_0").is_none() {
            p.phi0 = 0.863937979737193;
        }

        // if center long is not set use 42d30'E of Ferro - 17d40' for Ferro
        // that will correspond to using longitudes relative to greenwich
        // as input and output, instead of lat/long relative to Ferro
        if params.get("lon_0").is_none() {
            p.lam0 = 0.7417649320975901 - 0.308341501185665;
        }

        // if scale not set default to 0.9999
        if params.get("k").is_none() && params.get("k0").is_none() {
            p.k0 = 0.9999;
        }

        let easting_northing = !params.check_option("czech")?;

        // Set up shared parameters between forward and inverse
        let (e, es) = (p.ellps.e, p.ellps.es);
        let phi0 = p.phi0;
        let sinphi0 = phi0.sin();
        let alpha = (1. + (es * phi0.cos().powi(4)) / (1. - es)).sqrt();

        let u0 = (sinphi0 / alpha).asin();
        let g = ((1. + e * sinphi0) / (1. - e * sinphi0)).powf(alpha * e / 2.);

        let tan_half_phi0_plus_pi_4 = (phi0 / 2. + FRAC_PI_4).tan();
        if tan_half_phi0_plus_pi_4 == 0.0 {
            return Err(Error::InputStringError(
                "Invalid value for lat_0: lat_0 + PI/4 should be different from 0",
            ));
        }

        let n0 = (1. - es).sqrt() / (1. - es * sinphi0.powf(2.));

        Ok(Projection {
            e,
            xyfact: (2. * p.x0 / p.ellps.a, 2. * p.y0 / p.ellps.a),
            alpha,
            k: (u0 / 2. + FRAC_PI_4).tan() / tan_half_phi0_plus_pi_4.powf(alpha) * g,
            n: S0.sin(),
            rho0: p.k0 * n0 / S0.tan(),
            ad: FRAC_PI_2 - UQ,
            easting_northing,
        })
    }

    pub fn forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let sinphi = phi.sin();
        let gfi = ((1. + self.e * sinphi) / (1. - self.e * sinphi)).powf(self.alpha * self.e / 2.);
        let u = 2.
            * ((self.k * (phi / 2. + FRAC_PI_4).tan().powf(self.alpha) / gfi).atan() - FRAC_PI_4);

        let deltav = -lam * self.alpha;
        let s = (self.ad.cos() * u.sin() + self.ad.sin() * u.cos() * deltav.cos()).asin();
        let cos_s = s.cos();

        Ok(if cos_s < 1.0e-12 {
            (0., 0., z)
        } else {
            let eps = self.n * (u.cos() * deltav.sin() / cos_s).asin();
            let rho = self.rho0 * (S0 / 2. + FRAC_PI_4).tan().powf(self.n)
                / (s / 2. + FRAC_PI_4).tan().powf(self.n);

            let (x, y) = (rho * eps.sin(), rho * eps.cos());
            if self.easting_northing {
                (
                    // The default non-Czech convention uses easting, northing, so we have
                    // to reverse the sign of the coordinates. But to do so, we have to
                    // take into account the false easting/northing
                    -x - self.xyfact.0,
                    -y - self.xyfact.1,
                    z,
                )
            } else {
                (x, y, z)
            }
        })
    }

    fn inverse(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (x, y) = if self.easting_northing {
            // NOTE that correction factors are reversed in y/x
            (-y - self.xyfact.0, -x - self.xyfact.1)
        } else {
            (y, x)
        };

        let rho = x.hypot(y);
        let eps = y.atan2(x);

        let d = eps / S0.sin();
        let s = if rho == 0.0 {
            FRAC_PI_2
        } else {
            2. * (((self.rho0 / rho).powf(1. / self.n) * (S0 / 2. + FRAC_PI_4).tan()).atan()
                - FRAC_PI_4)
        };

        let u = (self.ad.cos() * s.sin() - self.ad.sin() * s.cos() * d.cos()).asin();
        let deltav = (s.cos() * d.sin() / u.cos()).asin();

        let lam = -deltav / self.alpha;

        let mut fi1 = u;
        let mut phi;
        for _ in 0..MAX_ITER {
            phi = 2.
                * ((self.k.powf(-1. / self.alpha)
                    * (u / 2. + FRAC_PI_4).tan().powf(1. / self.alpha)
                    * ((1. + self.e * fi1.sin()) / (1. - self.e * fi1.sin())).powf(self.e / 2.))
                .atan()
                    - FRAC_PI_4);
            if (fi1 - phi).abs() < EPS {
                return Ok((lam, phi, z));
            }
            fi1 = phi;
        }
        Err(Error::CoordTransOutsideProjectionDomain)
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

    // NOTE Krovak projection is valid for restricted bounding box
    // see https://epsg.io/5513

    #[test]
    fn proj_krovak() {
        let p = Proj::from_proj_string("+proj=krovak +units=m").unwrap();

        println!("{:#?}", p.projection());

        let inputs = [
            (
                (12.09, 47.73, 0.),
                (-951555.937880165293, -1276319.151569747366, 0.),
            ),
            (
                (22.56, 51.06, 0.),
                (-159523.534749580635, -983087.548008236452, 0.),
            ),
        ];

        test_proj_forward(&p, &inputs, 1e-6);
        test_proj_inverse(&p, &inputs, 1e-6);
    }
}
