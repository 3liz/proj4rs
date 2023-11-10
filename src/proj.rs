//!
//! Projection installation
//!
//! See <https://proj.org/usage/projections.html#cartographic-projection>
//! for parameter's descriptions.
//!

use crate::datum_params::DatumParams;
use crate::datum_transform::Datum;
use crate::datums::{self, DatumDefn};
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::parameters::ParamList;
use crate::projections::{find_projection, ProjDelegate};
use crate::{ellipsoids, prime_meridians, projstring, units};

use std::fmt;

pub type Axis = [u8; 3];

const NORMALIZED_AXIS: Axis = [b'e', b'n', b'u'];

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ProjType {
    Geocentric,
    Latlong,
    Other,
}

/// A Proj object hold informations and parameters
/// for a projection
#[derive(Debug, Clone)]
pub(crate) struct ProjData {
    pub(crate) ellps: Ellipsoid,
    pub(crate) axis: Axis,
    pub(crate) proj_type: ProjType,
    pub(crate) from_greenwich: f64, // prime meridian
    pub(crate) to_meter: f64,
    pub(crate) vto_meter: f64,
    pub(crate) x0: f64,
    pub(crate) y0: f64,
    pub(crate) k0: f64,
    pub(crate) lam0: f64,
    pub(crate) phi0: f64,
}

///
/// Projection definition
///
/// A projection is initalized with a projstring which is basically the same
/// as those used with proj4.
///
/// See <https://proj.org/usage/projections.html#cartographic-projection> for
/// a detailed description of the parameters
///
/// ```rust
/// use proj4rs::Proj;
///
/// // Create utm projection
/// let utm = Proj::from_proj_string("+proj=utm +ellps=GRS80 +zone=30").unwrap();
///
/// // Create latlon stub projection with ellipsoid "GRS80"
/// let geo = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
/// ```
#[derive(Clone)]
pub struct Proj {
    datum: Datum,
    geoc: bool,
    over: bool, // over-ranging flag
    // Units
    units: &'static str,
    vunits: &'static str,
    // Set by projections initialization
    projdata: ProjData,
    projname: &'static str,
    projection: ProjDelegate,
}

//----------------------
// Projection parameters
//----------------------
impl Proj {
    /// Return the projection name
    #[inline]
    pub fn projname(&self) -> &'static str {
        self.projname
    }
    #[inline]
    pub(crate) fn projection(&self) -> &ProjDelegate {
        &self.projection
    }
    /// Check if inverse projection exists
    #[inline]
    pub fn has_inverse(&self) -> bool {
        self.projection.has_inverse()
    }
    /// Check if forward projection exists
    #[inline]
    pub fn has_forward(&self) -> bool {
        self.projection.has_forward()
    }
    #[inline]
    pub(crate) fn data(&self) -> &ProjData {
        &self.projdata
    }
    #[inline]
    pub(crate) fn datum(&self) -> &Datum {
        &self.datum
    }
    #[inline]
    pub(crate) fn geoc(&self) -> bool {
        self.geoc
    }
    #[inline]
    pub(crate) fn over(&self) -> bool {
        self.over
    }
    // Delegate
    #[inline]
    pub(crate) fn ellipsoid(&self) -> &Ellipsoid {
        &self.projdata.ellps
    }
    #[inline]
    pub fn vto_meter(&self) -> f64 {
        self.projdata.vto_meter
    }
    #[inline]
    pub fn to_meter(&self) -> f64 {
        self.projdata.to_meter
    }
    #[inline]
    pub fn axis(&self) -> &Axis {
        &self.projdata.axis
    }
    /// Return true if the axis are normalized
    #[inline]
    pub fn is_normalized_axis(&self) -> bool {
        self.projdata.axis == NORMALIZED_AXIS
    }
    #[inline]
    pub fn is_latlong(&self) -> bool {
        self.projdata.proj_type == ProjType::Latlong
    }
    #[inline]
    pub fn is_geocent(&self) -> bool {
        self.projdata.proj_type == ProjType::Geocentric
    }
    #[inline]
    pub fn from_greenwich(&self) -> f64 {
        self.projdata.from_greenwich
    }

    #[inline]
    pub fn projection_type(&self) -> ProjType {
        self.projdata.proj_type
    }

    pub fn units(&self) -> &'static str {
        if self.is_latlong() {
            units::DEGREES
        } else {
            self.units
        }
    }

    #[inline]
    pub fn vunits(&self) -> &'static str {
        self.vunits
    }
}

//-------------------------
// Initialisation
//------------------------
impl Proj {
    // ----------------
    // Datum definition
    // ----------------
    fn get_datum_defn<'a>(params: &'a ParamList) -> Result<Option<&'a DatumDefn>> {
        // Do we have a "datum" parameter ?
        params
            .get("datum")
            .map(|p| match datums::find_datum(p.try_into()?) {
                Some(v) => Ok(Some(v)),
                None => Err(Error::InvalidDatum),
            })
            .unwrap_or(Ok(None))
    }

    // --------------
    // Prime meridian
    // --------------
    fn get_prime_meridian(params: &ParamList) -> Result<f64> {
        params
            .get("pm")
            .map(
                |p| match prime_meridians::find_prime_meridian(p.try_into()?) {
                    Some(v) => Ok(v),
                    None => p.try_into(),
                },
            )
            .unwrap_or(Ok(0.))
    }

    // -----------------
    // Datum parameters
    // ----------------
    fn get_datum_params(params: &ParamList, defn: Option<&DatumDefn>) -> Result<DatumParams> {
        // Precedence order is 'nadgrids', 'towgs84', 'datum'
        if let Some(p) = params.get("nadgrids") {
            // Nadgrids
            DatumParams::from_nadgrid_str(p.try_into()?)
        } else if let Some(p) = params.get("towgs84") {
            DatumParams::from_towgs84_str(p.try_into()?)
            // ToWGS84
        } else if let Some(p) = defn {
            DatumParams::try_from(&p.params)
        } else {
            Ok(DatumParams::default())
        }
    }

    // -----------------
    // Ellipsoid
    // ----------------
    fn get_ellipsoid(params: &ParamList, datum_def: Option<&DatumDefn>) -> Result<Ellipsoid> {
        if let Some(radius) = params.get("R") {
            // Sphere override everything
            Ellipsoid::sphere(radius.try_into()?)
        } else if let Some(p) = params.get("ellps") {
            // Return from ellipse definition
            match ellipsoids::find_ellipsoid(p.try_into()?) {
                Some(defn) => Ellipsoid::try_from_ellipsoid_with_params(defn, params),
                None => Err(Error::InvalidEllipsoid),
            }
        } else if let Some(defn) = datum_def {
            // Retrieve from datum definition + parameters
            Ellipsoid::try_from_ellipsoid_with_params(defn.ellps, params)
        } else if let Some(a) = params.get("a") {
            Ellipsoid::try_from_semi_major_axis(a.try_into()?, params)
        } else {
            // Get a free WGS84
            Ellipsoid::try_from_ellipsoid_with_params(&ellipsoids::constants::WGS84, params)
        }
    }

    // -----------------
    // Axis
    // ----------------
    fn get_axis(params: &ParamList) -> Result<Axis> {
        if let Some(p) = params.get("axis") {
            let axis_arg: &str = p.try_into()?;
            if axis_arg.len() != 3 {
                Err(Error::InvalidAxis)
            } else {
                let mut axis = [0u8, 0u8, 0u8];
                // Find Easting/Westing
                // This ensure that no token is repeated unless
                // one of the `find` will fail.
                let ew = axis_arg.find(['e', 'w']).ok_or(Error::InvalidAxis)?;
                let ns = axis_arg.find(['n', 's']).ok_or(Error::InvalidAxis)?;
                let ud = axis_arg.find(['u', 'd']).ok_or(Error::InvalidAxis)?;
                axis[ew] = axis_arg.as_bytes()[ew];
                axis[ns] = axis_arg.as_bytes()[ns];
                axis[ud] = axis_arg.as_bytes()[ud];
                Ok(axis)
            }
        } else {
            Ok(NORMALIZED_AXIS)
        }
    }

    // -----------------
    // Units
    // ----------------
    fn get_horizontal_units(params: &ParamList) -> Result<units::UnitDefn> {
        if let Some(p) = params.get("units") {
            let name: &str = p.try_into()?;
            if name.eq_ignore_ascii_case(units::DEGREES) {
                // Just a dummy value
                Ok(units::METER)
            } else {
                units::find_units(name).ok_or(Error::InvalidParameterValue("Invalid units"))
            }
        } else {
            Ok(params
                .try_value::<f64>("to_meter")?
                .map(units::from_value)
                .unwrap_or(units::METER))
        }
    }

    fn get_vertical_units(params: &ParamList) -> Result<units::UnitDefn> {
        if let Some(p) = params.get("vunits") {
            units::find_units(p.try_into()?).ok_or(Error::InvalidParameterValue("Invalid units"))
        } else {
            // XXX in proj4 vto_meter accept fractional expression: '/'
            Ok(params
                .try_value::<f64>("vto_meter")?
                .map(units::from_value)
                .unwrap_or(units::METER))
        }
    }

    ///
    /// Proj object constructor
    ///
    /// Consume a ParamList and create a Proj object
    ///
    pub fn init(params: ParamList) -> Result<Self> {
        // Find projection
        let proj_init = params
            .get("proj")
            .ok_or(Error::MissingProjectionError)
            .and_then(|name| find_projection(name.try_into()?).ok_or(Error::ProjectionNotFound))?;

        // Get datum definition (if any)
        let datum_defn = Self::get_datum_defn(&params)?;

        // Get datum parameters
        let datum_params = Self::get_datum_params(&params, datum_defn)?;

        // Do we have an ellipse ?
        let ellps = Self::get_ellipsoid(&params, datum_defn)?;

        // Get prime meridian
        let from_greenwich = Self::get_prime_meridian(&params)?;

        // Axis
        let axis = Self::get_axis(&params)?;

        // horizontal units
        let horz_units = Self::get_horizontal_units(&params)?;
        let vert_units = Self::get_vertical_units(&params)?;

        let to_meter = horz_units.to_meter;
        let vto_meter = vert_units.to_meter;

        // Datum
        let datum = Datum::new(&ellps, datum_params);

        let mut projdata = ProjData {
            ellps,
            axis,
            proj_type: ProjType::Other,
            from_greenwich,
            to_meter,
            vto_meter,
            // Central meridian_
            lam0: params.try_angular_value("lon_0")?.unwrap_or(0.),
            phi0: params.try_angular_value("lat_0")?.unwrap_or(0.),
            x0: params.try_value("x_0")?.unwrap_or(0.),
            y0: params.try_value("y_0")?.unwrap_or(0.),
            // Proj4 compatibility
            k0: match params.get("k0") {
                Some(p) => Some(p.try_into()).transpose(),
                None => params.try_value("k"),
            }?
            .unwrap_or(1.),
        };

        let project = proj_init.init(&mut projdata, &params)?;
        Ok(Self {
            datum,
            // Use Geocentric Latitude
            // see https://proj.org/operations/conversions/geoc.html
            geoc: params.check_option("geoc")?,
            over: params.check_option("over")?,
            units: horz_units.name,
            vunits: vert_units.name,
            projdata,
            projname: proj_init.name(),
            projection: project,
        })
    }

    /// Create from projstring definition
    pub fn from_proj_string(s: &str) -> Result<Self> {
        Self::init(projstring::parse(s)?)
    }

    /// Create projection from user string
    pub fn from_user_string(s: &str) -> Result<Self> {
        let s = s.trim();
        if s.starts_with('+') {
            Self::from_proj_string(s)
        } else if s.eq_ignore_ascii_case("WGS84") {
            Self::from_proj_string("+proj=longlat +ellps=WGS84")
        } else {
            Err(Error::UnrecognizedFormat)
        }
    }

    /// Create projection from user string
    ///
    /// # Examples
    ///
    /// ```
    /// use proj4rs::Proj;
    ///
    /// let proj = Proj::from_epsg_code(4326).unwrap();
    ///
    /// assert_eq!(
    ///     proj.projname(),
    ///     "longlat",
    /// );
    /// ```
    #[cfg(feature = "crs-definitions")]
    pub fn from_epsg_code(code: u16) -> Result<Self> {
        crs_definitions::from_code(code)
            .ok_or(Error::ProjectionNotFound)
            .and_then(|def| Self::from_proj_string(def.proj4))
    }
}

// -------------
// Display
// -------------
impl fmt::Debug for Proj {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "datum:      {:#?}", self.datum)?;
        writeln!(f, "geoc:       {:#?}", self.geoc)?;
        writeln!(f, "over:       {:#?}", self.over)?;
        writeln!(f, "data:       {:#?}", self.projdata)?;
        writeln!(f, "projname:   {:#?}", self.projname)?;
        writeln!(f, "projection: {:#?}", self.projection)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use crate::errors::{Error, Result};

    const INVALID_ELLPS: &str = "+proj=latlong +lon_0=5.937 +lat_ts=45.027 +ellps=foo";

    #[test]
    fn proj_invalid_ellps_param() {
        let p: Result<Proj> = Proj::from_proj_string(INVALID_ELLPS);

        assert!(p.is_err());
        let err = p.unwrap_err();
        println!("{:?}", err);
        assert!(matches!(err, Error::InvalidEllipsoid));
    }
}
