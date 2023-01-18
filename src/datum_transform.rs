//!
//! Datum transformation
//!
//! As with proj4/5 the datum transformation use WGS84 as hub for
//! converting data from one crs to another
//!
//! Datum shifts are carried out with the following steps:
//!  
//!    1. Convert (latitude, longitude, ellipsoidal height) to
//!       3D geocentric cartesian coordinates (X, Y, Z)
//!    2. Transform the (X, Y, Z) coordinates to the new datum, using a
//!       7 parameter Helmert transformation.
//!    3. Convert (X, Y, Z) back to (latitude, longitude, ellipsoidal height)
//!
//! Actually, the step 2 use WGS84 as conversion *hub* wich leads to apply
//! 2 Helmert transformations.
//!
//! With natgrids the steps are sligtly differents:
//!    1. Apply nadgrid transformation with source datum
//!    2. Convert to geocentric with source ellipsoid parameters
//!    3. Convert to geodetic with dest ellipsoid.
//!    4. Apply inverse nadgrids transformation vith destination datum
//!
use crate::datum_params::DatumParams;
use crate::ellps::Ellipsoid;
use crate::errors::Result;
use crate::nadgrids::NadgridShift;

pub mod utils {
    //!
    //! Geodetic calculation utilities
    //!
    use crate::constants::{FRAC_PI_2, PI, TAU};
    use crate::errors::{Error, Result};

    const GENAU: f64 = 1.0e-12;
    const GENAU2: f64 = GENAU * GENAU;
    const MAXITER: usize = 30;
    const FRAC_PI_2_EPS: f64 = 1.001 * FRAC_PI_2;

    /// Convert geodetic coordinates to geocentric coordinatesa
    ///
    /// The function Convert_Geodetic_To_Geocentric converts geodetic coordinates
    /// (latitude, longitude, and height) to geocentric coordinates (X, Y, Z),
    /// according to the current ellipsoid parameters.
    ///
    /// Latitude  : Geodetic latitude in radians                     (input)
    /// Longitude : Geodetic longitude in radians                    (input)
    /// Height    : Geodetic height, in meters                       (input)
    /// X         : Calculated Geocentric X coordinate, in meters    (output)
    /// Y         : Calculated Geocentric Y coordinate, in meters    (output)
    /// Z         : Calculated Geocentric Z coordinate, in meters    (output)
    pub fn geodetic_to_geocentric(
        x: f64,
        y: f64,
        z: f64,
        a: f64,
        es: f64,
    ) -> Result<(f64, f64, f64)> {
        let mut lon = x;
        let mut lat = y;

        if lat < -FRAC_PI_2 && lat > -FRAC_PI_2_EPS {
            lat = -FRAC_PI_2
        } else if lat > FRAC_PI_2 && lat < FRAC_PI_2_EPS {
            lat = FRAC_PI_2
        } else if !(-FRAC_PI_2..=FRAC_PI_2).contains(&lat) {
            return Err(Error::LatitudeOutOfRange);
        };

        if lon > PI {
            // TAU is 2PI
            lon -= TAU;
        }

        let sin_lat = lat.sin();
        // Earth radius at location
        let rn = a / (1. - es * (sin_lat * sin_lat)).sqrt();

        let cos_lat = lat.cos();
        Ok((
            (rn + z) * cos_lat * lon.cos(),
            (rn + z) * cos_lat * lon.sin(),
            ((rn * (1. - es)) + z) * sin_lat,
        ))
    }

    /// Convert geocentric coordinates to geodetic coordinates
    ///
    ///  ### Reference...
    ///
    /// Wenzel, H.-G.(1985): Hochauflösende Kugelfunktionsmodelle für
    /// das Gravitationspotential der Erde. Wiss. Arb. Univ. Hannover
    /// Nr. 137, p. 130-131.
    ///
    /// Programmed by GGA- Leibniz-Institute of Applied Geophysics
    ///               Stilleweg 2
    ///               D-30655 Hannover
    ///              Federal Republic of Germany
    ///              Internet: www.gga-hannover.de
    ///
    /// Hannover, March 1999, April 2004.
    /// see also: comments in statements
    ///
    /// remarks:
    /// Mathematically exact and because of symmetry of rotation-ellipsoid,
    /// each point (X,Y,Z) has at least two solutions (Latitude1,Longitude1,Height1) and
    /// (Latitude2,Longitude2,Height2). Is point=(0.,0.,Z) (P=0.), so you get even
    /// four solutions,»  every two symmetrical to the semi-minor axis.
    /// Here Height1 and Height2 have at least a difference in order of
    /// radius of curvature (e.g. (0,0,b)=> (90.,0.,0.) or (-90.,0.,-2b);
    /// (a+100.)*(sqrt(2.)/2.,sqrt(2.)/2.,0.) => (0.,45.,100.) or
    /// (0.,225.,-(2a+100.))).
    /// The algorithm always computes (Latitude,Longitude) with smallest |Height|.
    /// For normal computations, that means |Height|<10000.m, algorithm normally
    /// converges after to 2-3 steps!!!
    /// But if |Height| has the amount of length of ellipsoid's axis
    /// (e.g. -6300000.m),»   algorithm needs about 15 steps.
    pub fn geocentric_to_geodetic(
        x: f64,
        y: f64,
        z: f64,
        a: f64,
        es: f64,
        b: f64,
    ) -> Result<(f64, f64, f64)> {
        let d2 = (x * x) + (y * y);

        // distance between semi-minor axis and location
        let p = d2.sqrt();
        // distance between center and location
        let rr = (d2 + z * z).sqrt();

        // if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
        // of ellipsoid (=center of mass), Latitude becomes PI/2
        let lon = if p / a < GENAU {
            if rr / a < GENAU {
                return Ok((0., FRAC_PI_2, -b));
            }
            0.
        } else {
            y.atan2(x)
        };

        //--------------------------------------------------------------
        // Following iterative algorithm was developped by
        // Institut for Erdmessung", University of Hannover, July 1988.
        // Internet: www.ife.uni-hannover.de
        // Iterative computation of CPHI,SPHI and Height.
        // Iteration of CPHI and SPHI to 10**-12 radian resp.
        // 2*10**-7 arcsec.
        // --------------------------------------------------------------
        let ct = z / rr;
        let st = p / rr;
        let mut rx = 1.0 / (1.0 - es * (2.0 - es) * st * st).sqrt();
        let mut cphi0 = st * (1.0 - es) * rx;
        let mut sphi0 = ct * rx;
        let (mut rk, mut rn, mut cphi, mut sphi, mut sdphi, mut height);

        // loop to find sin(Latitude) resp. Latitude
        // until |sin(Latitude(iter)-Latitude(iter-1))| < genau

        // Note: using `for _ in 0..MAXITER { ... }` lead to compiler error
        // about unitialized variables
        let mut iter = 0;
        loop {
            iter += 1;
            rn = a / (1.0 - es * sphi0 * sphi0).sqrt();
            // ellipsoidal (geodetic) height
            height = p * cphi0 + z * sphi0 - rn * (1.0 - es * sphi0 * sphi0);

            //  avoid zero division
            if (rn + height) == 0. {
                return Ok((lon, 0., height));
            }

            rk = es * rn / (rn + height);
            rx = 1.0 / (1.0 - rk * (2.0 - rk) * st * st).sqrt();
            cphi = st * (1.0 - rk) * rx;
            sphi = ct * rx;
            sdphi = sphi * cphi0 - cphi * sphi0;
            cphi0 = cphi;
            sphi0 = sphi;

            if sdphi * sdphi <= GENAU2 {
                break;
            }

            if iter >= MAXITER {
                break;
            }
        }

        // ellipsoidal (geodetic) latitude
        Ok((lon, sphi.atan2(cphi.abs()), height))
    }
}

use utils::*;
use DatumParams::*;

/// Hold datum Informations
#[derive(Default, Debug)]
pub struct Datum<GS: NadgridShift> {
    params: DatumParams<GS>,
    a: f64,
    b: f64,
    es: f64,
}

impl<GS: NadgridShift> Datum<GS> {
    pub fn new(ellps: &Ellipsoid, params: DatumParams<GS>) -> Self {
        Self {
            // check for WGS84/GRS80
            params: if params == ToWGS84_3(0., 0., 0.)
                && ellps.a == 6378137.0
                && (ellps.es - 0.006694379990).abs() < 0.000000000050
            {
                ToWGS84_0
            } else {
                params
            },
            a: ellps.a,
            b: ellps.b,
            es: ellps.es,
        }
    }

    /// Convert from geodetic coordinates to wgs84/geocentric
    #[inline(always)]
    fn towgs84(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match &self.params {
            ToWGS84_0 => geodetic_to_geocentric(x, y, z, self.a, self.es),
            ToWGS84_3(dx, dy, dz) => geodetic_to_geocentric(x, y, z, self.a, self.es)
                .map(|(x, y, z)| (x + dx, y + dy, z + dz)),
            ToWGS84_7(dx, dy, dz, rx, ry, rz, s) => {
                geodetic_to_geocentric(x, y, z, self.a, self.es).map(|(x, y, z)| {
                    (
                        dx + s * (x - rz * y + ry * z),
                        dy + s * (rz * x + y - rx * z),
                        dz + s * (-ry * x + rx * y + z),
                    )
                })
            }
            NadGrids(grids) => grids
                .apply(false, x, y, z)
                .and_then(|(x, y, z)| geodetic_to_geocentric(x, y, z, self.a, self.es)),
            NoDatum => Ok((x, y, z)),
        }
    }

    /// Convert from geocentric/wgs84 to geodetic coordinates
    #[inline(always)]
    fn fromwgs84(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        match &self.params {
            ToWGS84_0 => geocentric_to_geodetic(x, y, z, self.a, self.es, self.b),
            ToWGS84_3(dx, dy, dz) => {
                geocentric_to_geodetic(x - dx, y - dy, z - dz, self.a, self.es, self.b)
            }
            ToWGS84_7(dx, dy, dz, rx, ry, rz, s) => {
                let (x, y, z) = ((x - dx) / s, (x - dy) / s, (y - dz) / s);
                geocentric_to_geodetic(
                    x + rz * y - ry * z,
                    -rz * x + y + rx * z,
                    ry * x - rx * y + z,
                    self.a,
                    self.es,
                    self.b,
                )
            }
            NadGrids(grids) => geocentric_to_geodetic(x, y, y, self.a, self.es, self.b)
                .and_then(|(x, y, z)| grids.apply(true, x, y, z)),
            NoDatum => Ok((x, y, z)),
        }
    }

    pub fn use_nadgrids(&self) -> bool {
        matches!(self.params, NadGrids(_))
    }

    pub fn no_datum(&self) -> bool {
        matches!(self.params, NoDatum)
    }

    pub fn has_wgs84_params(&self) -> bool {
        matches!(self.params, ToWGS84_0 | ToWGS84_3(..) | ToWGS84_7(..))
    }

    pub fn is_identical_to(&self, other: &Self) -> bool {
        // the tolerance for es is to ensure that GRS80 and WGS84
        // are considered identical
        (self.use_nadgrids() && other.use_nadgrids() || self.params == other.params)
            && self.a == other.a
            && (self.es - other.es).abs() < 0.000000000050
    }
}

const SRS_WGS84_SEMIMAJOR: f64 = 6378137.0;
const SRS_WGS84_SEMIMINOR: f64 = 6356752.314;
const SRS_WGS84_ES: f64 = 0.0066943799901413165;

/// Datum transform
pub struct DatumTransform<G: NadgridShift> {
    pub src: Datum<G>,
    pub dst: Datum<G>,
    pub identity: bool,
}

impl<G: NadgridShift> DatumTransform<G> {
    pub fn new(mut src: Datum<G>, mut dst: Datum<G>) -> Self {
        // Change ellipse parameters to wgs84
        // when using nadgrids
        if src.use_nadgrids() {
            src.a = SRS_WGS84_SEMIMAJOR;
            src.b = SRS_WGS84_SEMIMINOR;
            src.es = SRS_WGS84_ES;
        }

        if dst.use_nadgrids() {
            dst.a = SRS_WGS84_SEMIMAJOR;
            dst.b = SRS_WGS84_SEMIMINOR;
            dst.es = SRS_WGS84_ES;
        }

        let identity = src.params == NoDatum
            || dst.params == NoDatum
            || src.is_identical_to(&dst)
            || src.a == dst.a
                && src.es == dst.es
                && !src.has_wgs84_params()
                && !dst.has_wgs84_params();

        Self { src, dst, identity }
    }

    #[inline(always)]
    pub fn transform(&self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        if self.identity {
            Ok((x, y, z))
        } else {
            self.src
                .towgs84(x, y, z)
                .and_then(|(x, y, z)| self.dst.fromwgs84(x, y, z))
        }
    }
}
