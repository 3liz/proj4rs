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
use crate::geocent::{geocentric_to_geodetic, geodetic_to_geocentric};
use crate::transform::Direction;

use DatumParams::*;

const SRS_WGS84_SEMIMAJOR: f64 = 6378137.0;
const SRS_WGS84_SEMIMINOR: f64 = 6356752.314;
const SRS_WGS84_ES: f64 = 0.0066943799901413165;

/// Hold datum Informations
#[derive(Debug)]
pub(crate) struct Datum {
    params: DatumParams,
    pub a: f64,
    pub b: f64,
    pub es: f64,
}

impl Datum {
    pub fn new(ellps: &Ellipsoid, params: DatumParams) -> Self {
        // Change ellipse parameters to wgs84
        // when using nadgrids
        let (a, b, es) = if params.use_nadgrids() {
            (SRS_WGS84_SEMIMAJOR, SRS_WGS84_SEMIMINOR, SRS_WGS84_ES)
        } else {
            (ellps.a, ellps.b, ellps.es)
        };

        Self {
            // check for WGS84/GRS80
            params: if params == ToWGS84_3(0., 0., 0.)
                && ellps.a == SRS_WGS84_SEMIMAJOR
                && (ellps.es - SRS_WGS84_ES).abs() < 0.000000000050
            {
                ToWGS84_0
            } else {
                params
            },
            a,
            b,
            es,
        }
    }

    /// Convert from geodetic coordinates to wgs84/geocentric
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
                .apply_shift(Direction::Forward, x, y, z)
                .and_then(|(x, y, z)| geodetic_to_geocentric(x, y, z, self.a, self.es)),
            NoDatum => Ok((x, y, z)),
        }
    }

    /// Convert from geocentric/wgs84 to geodetic coordinates
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
                .and_then(|(x, y, z)| grids.apply_shift(Direction::Inverse, x, y, z)),
            NoDatum => Ok((x, y, z)),
        }
    }

    #[inline]
    pub fn use_nadgrids(&self) -> bool {
        self.params.use_nadgrids()
    }

    #[inline]
    pub fn no_datum(&self) -> bool {
        self.params.no_datum()
    }

    #[inline]
    pub fn use_towgs84(&self) -> bool {
        self.params.use_towgs84()
    }

    /// Return true if the datum are identical
    pub fn is_identical_to(&self, other: &Self) -> bool {
        // the tolerance for es is to ensure that GRS80 and WGS84
        // are considered identical
        (self.params == other.params)
            && self.a == other.a
            && (self.es - other.es).abs() < 0.000000000050
    }

    /// Transform geographic coordinates between datums
    ///
    /// No identity checking is done
    #[inline]
    pub fn transform(src: &Self, dst: &Self, x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
        src.towgs84(x, y, z)
            .and_then(|(x, y, z)| dst.fromwgs84(x, y, z))
    }
}
