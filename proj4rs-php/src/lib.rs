//!
//! Php binding entry point
//!
//! See https://davidcole1340.github.io/ext-php-rs


use proj4rs::{errors, proj, transform};
use ext_php_rs::prelude::*;


#[cfg(feature = "logging")]
use log;

// Entry point
//pub fn main() {
//    #[cfg(feature = "logging")]
//    console_log::init_with_level(log::Level::Trace).unwrap();
//}

// ----------------------------
// Wrapper for Projection
// ---------------------------
#[php_class]
pub struct Projection {
    inner: proj::Proj,
}

impl From<proj::Proj> for Projection {
    fn from(p: proj::Proj) -> Self {
        Self { inner: p }
    }
}


#[php_impl(rename_methods = "camelCase")]
impl Projection {

    #[constructor]
    fn new(defn: &str) -> PhpResult<Self> {
     proj::Proj::from_user_string(defn)
        .map(Projection::from)
        .map_err(|e| PhpException::from(e.to_string()))
    }

    // see https://github.com/davidcole1340/ext-php-rs/issues/325   
    // pub fn projname(&self) -> &'static str {
    //   self.inner.projname()
    // }
    #[getter(rename = "projName")]
    pub fn projname(&self) -> String {
        self.inner.projname().into()
    }

    #[getter(rename = "isLatlong")]
    pub fn is_latlong(&self) -> bool {
        self.inner.is_latlong()
    }

    #[getter(rename = "isGeocentric")]
    pub fn is_geocent(&self) -> bool {
        self.inner.is_geocent()
    }

    #[getter]
    pub fn axis(&self) -> String {
        String::from_utf8_lossy(self.inner.axis()).into_owned()
    }

    #[getter(rename = "isNormalizedAxis")]
    pub fn is_normalized_axis(&self) -> bool {
        self.inner.is_normalized_axis()
    }

    #[getter(rename = "toMeter")]
    pub fn to_meter(&self) -> f64 {
        self.inner.to_meter()
    }

    // see https://github.com/davidcole1340/ext-php-rs/issues/325   
    // pub fn units(&self) -> &'static str {
    //    self.inner.units()
    // }
    #[getter]
    pub fn units(&self) -> String {
        self.inner.units().into()
    }
}


// ----------------------------
// Wrapper for Transform
// ---------------------------
#[php_class]
pub struct Point {
    #[prop]
    pub x: f64,
    #[prop]
    pub y: f64,
    #[prop]
    pub z: f64,
}

#[php_impl]
impl Point {
    pub fn __construct(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

// Transform Point

impl transform::Transform for Point {
    /// Strict mode: return exception
    /// as soon as with have invalid coordinates or
    /// that the reprojection failed
    fn transform_coordinates<F>(&mut self, f: &mut F) -> errors::Result<()>
    where
        F: FnMut(f64, f64, f64) -> errors::Result<(f64, f64, f64)>,
    {
        f(self.x, self.y, self.z).map(|(x, y, z)| {
            self.x = x;
            self.y = y;
            self.z = z;
        })
    }
}

#[php_function]
pub fn transform_point(
    src: &Projection,
    dst: &Projection,
    point: &mut Point,
    convert: bool,
) -> PhpResult<()> {
    if point.x.is_nan() || point.y.is_nan() {
        return Err(PhpException::from(errors::Error::NanCoordinateValue.to_string()));
    }

    if convert && src.inner.is_latlong() {
        point.x = point.x.to_radians();
        point.y = point.y.to_radians();
    }

    transform::transform(&src.inner, &dst.inner, point)
        .map_err(|e| PhpException::from(e.to_string()))?;

    if convert && dst.inner.is_latlong() {
        point.x = point.x.to_degrees();
        point.y = point.y.to_degrees();
    }
    Ok(())
}


#[php_module]
pub fn get_module(module: ModuleBuilder) -> ModuleBuilder {
    module
}
