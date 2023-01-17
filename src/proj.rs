//!
//! Projection installation
//!

use crate::datum_params::DatumParams;
use crate::datum_transform::Datum;
use crate::datums::{self, DatumDefn};
use crate::ellps::Ellipsoid;
use crate::errors::{Error, Result};
use crate::nadgrids::{NadgridShift, NullGridShift};
use crate::parameters::ParamList;
use crate::{ellipsoids, prime_meridians, projstring};

//============================
//
// Projection
//
//============================

pub struct Projection<GS: NadgridShift = NullGridShift> {
    pj_pm: f64,
    pj_ellps: Ellipsoid,
    pj_datum: Datum<GS>,
}

impl<GS: NadgridShift> Projection<GS> {
    // ----------------
    // Datum definition
    // ----------------
    fn datum_defn<'a>(params: &'a ParamList) -> Result<Option<&'a DatumDefn>> {
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
    fn prime_meridian(params: &ParamList) -> Result<f64> {
        params
            .get("pm")
            .map(
                |p| match prime_meridians::find_prime_meridian(p.try_into()?) {
                    Some(v) => Ok(v),
                    None => p.try_convert::<f64>(),
                },
            )
            .unwrap_or(Ok(0.))
    }

    // -----------------
    // Datum parameters
    // ----------------
    fn datum_params(params: &ParamList, defn: Option<&DatumDefn>) -> Result<DatumParams<GS>> {
        // Precedence order is 'nadgrids', 'towgs84', 'datum'
        if let Some(p) = params.get("nadgrids") {
            // Nadgrids
            DatumParams::from_nagrid_str(p.try_into()?)
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
    fn ellipsoid(params: &ParamList, datum_def: Option<&DatumDefn>) -> Result<Ellipsoid> {
        if let Some(radius) = params.get("R") {
            // Sphere override everything
            Ellipsoid::sphere(radius.try_into()?)
        } else if let Some(p) = params.get("ellps") {
            // Return from ellipse definition
            match ellipsoids::find_ellipsoid(p.try_into()?) {
                Some(defn) => Ellipsoid::try_from_ellipsoid(defn),
                None => Err(Error::InvalidEllipsoid),
            }
        } else if let Some(defn) = datum_def {
            // Retrieve from datum definition + parameters
            Ellipsoid::try_from_ellipsoid_with_params(defn.ellps, params)
        } else {
            // Get a free WGS84
            Ok(Ellipsoid::default())
        }
    }

    /// Consume a ParamList and create a Projection object
    pub fn init(params: ParamList) -> Result<Self> {
        // Projection name
        let _projname = params.get("proj").ok_or(Error::MissingProjectionError);

        // Get datum definition (if any)
        let datum_defn = Self::datum_defn(&params)?;

        // Get datum parameters
        let datum_params = Self::datum_params(&params, datum_defn)?;

        // Do we have an ellipse ?
        let pj_ellps = Self::ellipsoid(&params, datum_defn)?;

        // Get prime meridian
        let pj_pm = Self::prime_meridian(&params)?;

        // Datum
        let pj_datum = Datum::new(&pj_ellps, datum_params);

        Ok(Self {
            pj_pm,
            pj_ellps,
            pj_datum,
        })
    }

    /// Create projection from string
    pub fn from_projstr(s: &str) -> Result<Self> {
        Self::init(projstring::parse(s)?)
    }
}
