//!
//! Geodetic to/from geocentrique  conversion
//!
use crate::errors::{Error, Result};
use crate::math::consts::{FRAC_PI_2, PI, TAU};

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
///
/// This conversion converts geodetic coordinate values (longitude, latitude, elevation above ellipsoid)
/// to their geocentric (X, Y, Z) representation, where the first axis (X) points from the Earth centre
/// to the point of longitude=0, latitude=0, the second axis (Y) points from the
/// Earth centre to the point of longitude=90, latitude=0 and the third axis (Z) points to the North pole
///
pub fn geodetic_to_geocentric(x: f64, y: f64, z: f64, a: f64, es: f64) -> Result<(f64, f64, f64)> {
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

    let (sin_lat, cos_lat) = lat.sin_cos();
    // Earth radius at location
    let rn = a / (1. - es * (sin_lat * sin_lat)).sqrt();
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
    // Following iterative algorithm was developed by
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
