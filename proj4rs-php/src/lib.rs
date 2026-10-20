//!
//! PHP extension for proj4rs
//!
//! Exposes the same API as the FFI bindings in `proj4rs-clib/php`:
//! `Proj4rs\Proj`, `Proj4rs\Transform` and `Proj4rs\Proj4rsException`.
//!
#![cfg_attr(windows, feature(abi_vectorcall))]
use ext_php_rs::{
    convert::IntoZval,
    exception::PhpException,
    flags::ClassFlags,
    prelude::*,
    types::{ZendClassObject, ZendHashTable, Zval},
    zend::ClassEntry,
};

use proj4rs::{
    errors, proj,
    transform::{TransformClosure, transform},
};

// SPL exception classes, exported by the PHP binary
unsafe extern "C" {
    static spl_ce_RuntimeException: *mut ClassEntry;
    static spl_ce_InvalidArgumentException: *mut ClassEntry;
}

fn runtime_exception() -> &'static ClassEntry {
    unsafe { spl_ce_RuntimeException.as_ref() }.expect("RuntimeException class not initialized")
}

fn invalid_argument(msg: &str) -> PhpException {
    let ce = unsafe { spl_ce_InvalidArgumentException.as_ref() }
        .expect("InvalidArgumentException class not initialized");
    PhpException::new(msg.into(), 0, ce)
}

fn proj4rs_error(err: errors::Error) -> PhpException {
    PhpException::from_class::<Proj4rsException>(err.to_string())
}

#[php_class]
#[php(name = "Proj4rs\\Proj4rsException")]
#[php(extends(ce = runtime_exception, stub = "\\RuntimeException"))]
#[derive(Default)]
pub struct Proj4rsException;

// ----------------------------
// Proj
// ----------------------------

/// A projection definition
#[php_class]
#[php(name = "Proj4rs\\Proj")]
#[php(flags = ClassFlags::Final)]
#[derive(Clone)]
pub struct Proj(proj::Proj);

#[php_impl]
impl Proj {
    /// `defn` may be:
    ///
    /// - A projstring
    /// - "WGS84" - equivalent to the projstring "+proj=longlat +ellps=WGS84"
    /// - An EPSG code (i.e "EPSG:3035") if the extension is built
    ///   with the `crs-definitions` feature
    pub fn __construct(defn: &str) -> PhpResult<Self> {
        proj::Proj::from_user_string(defn)
            .map(Self)
            .map_err(proj4rs_error)
    }

    pub fn projname(&self) -> &'static str {
        self.0.projname()
    }

    pub fn is_latlong(&self) -> bool {
        self.0.is_latlong()
    }

    pub fn is_geocent(&self) -> bool {
        self.0.is_geocent()
    }

    /// Return the axis directions as a 3 chars string (i.e "enu")
    pub fn axis(&self) -> String {
        String::from_utf8_lossy(self.0.axis()).into_owned()
    }

    pub fn is_normalized_axis(&self) -> bool {
        self.0.is_normalized_axis()
    }

    pub fn to_meter(&self) -> f64 {
        self.0.to_meter()
    }

    /// Return units of the projection (i.e "degrees", "m", "km", ...)
    pub fn units(&self) -> &'static str {
        self.0.units()
    }
}

// ----------------------------
// Transform
// ----------------------------

/// Coordinate transformation between two projections
///
/// Coordinates that fail to be projected are set to NAN.
#[php_class]
#[php(name = "Proj4rs\\Transform")]
#[php(flags = ClassFlags::Final)]
pub struct Transform {
    src: Proj,
    dst: Proj,
}

/// Accept either a `Proj` object or a definition string
fn to_proj(arg: &Zval) -> PhpResult<Proj> {
    if let Some(defn) = arg.str() {
        Proj::__construct(defn)
    } else if let Some(obj) = arg.extract::<&ZendClassObject<Proj>>() {
        Ok((**obj).clone())
    } else {
        Err(invalid_argument(
            "Expecting a Proj4rs\\Proj object or a string",
        ))
    }
}

/// Points as a flat interleaved buffer of `dim` coordinates per point
struct Points<'a> {
    buf: &'a mut [f64],
    dim: usize,
}

impl proj4rs::transform::Transform for Points<'_> {
    fn transform_coordinates<F: TransformClosure>(&mut self, f: &mut F) -> errors::Result<()> {
        let dim = self.dim;
        self.buf.chunks_exact_mut(dim).for_each(|p| {
            let z = if dim == 3 { p[2] } else { 0. };
            match f(p[0], p[1], z) {
                Ok((x, y, z)) => {
                    p[0] = x;
                    p[1] = y;
                    if dim == 3 {
                        p[2] = z;
                    }
                }
                Err(_) => p.fill(f64::NAN),
            }
        });
        Ok(())
    }
}

fn to_f64(v: &Zval) -> PhpResult<f64> {
    v.double()
        .or_else(|| v.long().map(|l| l as f64))
        .ok_or_else(|| invalid_argument("Expecting numeric coordinates"))
}

fn to_vec(arr: &ZendHashTable) -> PhpResult<Vec<f64>> {
    arr.values().map(to_f64).collect()
}

impl Transform {
    /// Transform a flat interleaved buffer in place
    fn apply(&self, buf: &mut [f64], dim: usize, convert: bool) -> PhpResult<()> {
        let (src, dst) = (&self.src.0, &self.dst.0);
        if convert && src.is_latlong() {
            buf.chunks_exact_mut(dim).for_each(|p| {
                p[0] = p[0].to_radians();
                p[1] = p[1].to_radians();
            });
        }
        transform(
            src,
            dst,
            &mut Points {
                buf: &mut *buf,
                dim,
            },
        )
        .map_err(proj4rs_error)?;
        if convert && dst.is_latlong() {
            buf.chunks_exact_mut(dim).for_each(|p| {
                p[0] = p[0].to_degrees();
                p[1] = p[1].to_degrees();
            });
        }
        Ok(())
    }
}

fn into_zval<T: IntoZval>(v: T) -> PhpResult<Zval> {
    v.into_zval(false).map_err(|e| e.to_string().into())
}

#[php_impl]
impl Transform {
    pub fn __construct(src: &Zval, dst: &Zval) -> PhpResult<Self> {
        Ok(Self {
            src: to_proj(src)?,
            dst: to_proj(dst)?,
        })
    }

    pub fn source(&self) -> Proj {
        self.src.clone()
    }

    pub fn destination(&self) -> Proj {
        self.dst.clone()
    }

    /// Transform coordinates given as separate scalars or arrays of values.
    ///
    /// If `convert` is true, latlong coordinates are assumed to be in degrees.
    ///
    /// Returns `[x, y]` or `[x, y, z]`, as floats if the inputs are scalars,
    /// as lists of floats otherwise.
    #[php(defaults(convert = true))]
    pub fn transform(
        &self,
        x: &Zval,
        y: &Zval,
        z: Option<&Zval>,
        convert: bool,
    ) -> PhpResult<Vec<Zval>> {
        let z = z.filter(|z| !z.is_null());
        let coords: Vec<&Zval> = [Some(x), Some(y), z].into_iter().flatten().collect();
        let dim = coords.len();

        let scalar = x.array().is_none();
        if coords.iter().any(|c| c.array().is_none() != scalar) {
            return Err(invalid_argument(
                "Coordinates must be all scalars or all arrays",
            ));
        }

        let columns = coords
            .iter()
            .map(|c| match c.array() {
                Some(arr) => to_vec(arr),
                None => to_f64(c).map(|v| vec![v]),
            })
            .collect::<PhpResult<Vec<_>>>()?;

        let len = columns[0].len();
        if columns.iter().any(|c| c.len() != len) {
            return Err(invalid_argument(
                "Coordinate arrays must have the same length",
            ));
        }

        let mut buf: Vec<f64> = (0..len)
            .flat_map(|i| columns.iter().map(move |c| c[i]))
            .collect();
        self.apply(&mut buf, dim, convert)?;

        (0..dim)
            .map(|k| {
                if scalar {
                    into_zval(buf[k])
                } else {
                    into_zval(buf.iter().skip(k).step_by(dim).copied().collect::<Vec<_>>())
                }
            })
            .collect()
    }

    /// Transform a single point given as `[x, y]` or `[x, y, z]`
    #[php(defaults(convert = true))]
    pub fn transform_point(&self, point: &ZendHashTable, convert: bool) -> PhpResult<Vec<f64>> {
        let mut buf = to_vec(point)?;
        let dim = buf.len();
        if dim != 2 && dim != 3 {
            return Err(invalid_argument("Points must have 2 or 3 coordinates"));
        }
        self.apply(&mut buf, dim, convert)?;
        Ok(buf)
    }

    /// Transform a list of points, each given as `[x, y]` or `[x, y, z]`.
    ///
    /// All points must have the same dimension.
    #[php(defaults(convert = true))]
    pub fn transform_points(
        &self,
        points: &ZendHashTable,
        convert: bool,
    ) -> PhpResult<Vec<Vec<f64>>> {
        let mut dim = 0;
        let mut buf = Vec::with_capacity(points.len() * 3);
        for p in points.values() {
            let p = p
                .array()
                .ok_or_else(|| invalid_argument("Points must be arrays"))?;
            let n = p.len();
            if dim == 0 {
                if n != 2 && n != 3 {
                    return Err(invalid_argument("Points must have 2 or 3 coordinates"));
                }
                dim = n;
            } else if n != dim {
                return Err(invalid_argument("All points must have the same dimension"));
            }
            for v in p.values() {
                buf.push(to_f64(v)?);
            }
        }
        if dim == 0 {
            return Ok(vec![]);
        }
        self.apply(&mut buf, dim, convert)?;
        Ok(buf.chunks_exact(dim).map(<[f64]>::to_vec).collect())
    }
}

#[php_module]
pub fn get_module(module: ModuleBuilder) -> ModuleBuilder {
    module
        .class::<Proj4rsException>()
        .class::<Proj>()
        .class::<Transform>()
}
