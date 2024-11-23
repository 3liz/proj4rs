use libc::{c_char, c_int, ssize_t};
use std::cell::RefCell;
use std::ffi::{CStr, CString};
use std::ptr;

use proj4rs::{
    errors, proj,
    transform::{transform, Transform, TransformClosure},
};

thread_local!(static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("").unwrap()));

fn set_last_error(err: errors::Error) {
    let error = err.to_string();
    LAST_ERROR.with(|c| {
        c.replace(CString::new(error.as_str()).unwrap_or(CString::new("Unknown").unwrap()))
    });
}

fn to_c_unit(name: &str) -> &'static str {
    match name {
        "km" => "km\0",
        "m" => "m\0",
        "dm" => "dm\0",
        "cm" => "cm\0",
        "mm" => "mm\0",
        "kmi" => "kmi\0",
        "in" => "in\0",
        "ft" => "ft\0",
        "yd" => "yd\0",
        "mi" => "mi\0",
        "fath" => "fath\0",
        "ch" => "ch\0",
        "link" => "link\0",
        "us-in" => "us-in\0",
        "us-ft" => "us-ft\0",
        "us-yd" => "us-yd\0",
        "us-ch" => "us-ch\0",
        "us-mi" => "us-mi\0",
        "ind-yd" => "ind-yd\0",
        "ind-ft" => "ind-ft\0",
        "ind-ch" => "ind-ch\0",
        "degrees" => "degrees\0",
        _ => "\0",
    }
}

/// Return the last error message
#[no_mangle]
pub extern "C" fn proj4rs_last_error() -> *const c_char {
    LAST_ERROR.with_borrow(|c| c.as_ptr() as *const c_char)
}

/// Opaque structure holding the internal representation
/// of projection.
pub struct Proj4rs {
    inner: proj::Proj,
    name: CString,
}

/// Create new projection object
///
/// `c_defn` may be:
///
/// * A projstring
/// * A "WGS84" string - equivalent to the projstring "+proj=longlat +ellps=WGS84"
/// * An EPSG code
///
#[no_mangle]
pub extern "C" fn proj4rs_proj_new(c_defn: *const c_char) -> *mut Proj4rs {
    let cstr_defn = unsafe { CStr::from_ptr(c_defn) };
    match cstr_defn
        .to_str()
        .map_err(errors::Error::from)
        .and_then(proj::Proj::from_user_string)
        .map(|p| Proj4rs {
            name: CString::new(p.projname()).unwrap(),
            inner: p,
        }) {
        Ok(p) => Box::into_raw(Box::new(p)),
        Err(err) => {
            set_last_error(err);
            ptr::null_mut()
        }
    }
}

/// Delete projection object
#[no_mangle]
pub extern "C" fn proj4rs_proj_delete(c_ptr: *mut Proj4rs) {
    if !c_ptr.is_null() {
        unsafe {
            let _ = Box::from_raw(c_ptr);
        }
    }
}

/// Returns the projection name
#[no_mangle]
pub extern "C" fn proj4rs_proj_projname(c_ptr: *const Proj4rs) -> *const c_char {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.name.as_ptr() as *const c_char
}

/// Returns true if the projection is geographic
#[no_mangle]
pub extern "C" fn proj4rs_proj_is_latlong(c_ptr: *const Proj4rs) -> bool {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.inner.is_latlong()
}

/// Returns true if the projection is geocentric
#[no_mangle]
pub extern "C" fn proj4rs_proj_is_geocent(c_ptr: *const Proj4rs) -> bool {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.inner.is_geocent()
}

/// Return the projection axes
///
/// The value returned is a pointer to 3-value byte array
///
/// The first value represent the direction `x` axis: 'e' (East) or 'w' (West)
/// The second value represent the direction of the `y` axis: 'n' (North) or 's' (South)
/// The third value represent the direction of  the `z` axis: 'u' (Up) or 'd' (Down)
///
/// Example:  
/// ```C
/// axes = proj4rs_proj_axis(Proj4rs);
/// if( axes[0] == 'e' && axes[1] == 'n' && axes[2] == 'u') {
///    printf("Axis are normalized\n")
/// }
/// ```
///
#[no_mangle]
pub extern "C" fn proj4rs_proj_axis(c_ptr: *const Proj4rs) -> *const u8 {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.inner.axis().as_ptr()
}

/// Return true if the axis are noramilized
#[no_mangle]
pub extern "C" fn proj4rs_proj_is_normalized_axis(c_ptr: *const Proj4rs) -> bool {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.inner.is_normalized_axis()
}

#[no_mangle]
pub extern "C" fn proj4rs_proj_to_meter(c_ptr: *const Proj4rs) -> f64 {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    proj.inner.to_meter()
}

/// Return units of the projection (i.e "degrees", "m", "km", ...)
#[no_mangle]
pub extern "C" fn proj4rs_proj_units(c_ptr: *const Proj4rs) -> *const c_char {
    assert!(!c_ptr.is_null(), "Null proj pointer");
    let proj: &Proj4rs = unsafe { &*c_ptr };
    to_c_unit(proj.inner.units()).as_ptr() as *const c_char
}

// ----------------------------
// Wrapper for Transform
// ---------------------------

pub const OK: c_int = 1;
pub const ERR: c_int = 0;

/// Transform array of coordinates from src
///
/// Projected coordinates will be set to  NAN in case of projection failure
///
/// * `x`, `y`, `z` are pointers to double.
/// * `z` may be null
/// * `x` and `y` and `z' - if not null - must hold the same number of values
/// * `stride` is the memory offset between two consecutive values
/// * `len` in the number of points
///
/// If `convert` is `true` then latlong coordinates are assumed te be in degrees.
///
#[no_mangle]
pub extern "C" fn proj4rs_transform(
    src_ptr: *const Proj4rs,
    dst_ptr: *const Proj4rs,
    x: *mut f64,
    y: *mut f64,
    z: *mut f64,
    len: ssize_t,
    stride: ssize_t,
    convert: bool,
) -> c_int {
    if src_ptr.is_null()
        || dst_ptr.is_null()
        || x.is_null()
        || y.is_null()
        || len <= 0
        || stride <= 0
    {
        return ERR;
    }

    let src: &Proj4rs = unsafe { &*src_ptr };
    let dst: &Proj4rs = unsafe { &*dst_ptr };

    if convert && src.inner.is_latlong() {
        unsafe {
            to_radians(x, y, len, stride);
        }
    }
    
    if let Err(err) = transform(&src.inner, &dst.inner, &mut Coords(x, y, z, len, stride)) {
        set_last_error(err);
        ERR
    } else {
        if convert && dst.inner.is_latlong() {
            unsafe {
                to_degrees(x, y, len, stride);
            }
        }
        OK
    }
}

struct Coords(*mut f64, *mut f64, *mut f64, isize, isize);

impl Transform for Coords {
    fn transform_coordinates<F: TransformClosure>(&mut self, f: &mut F) -> errors::Result<()>
    where
        F: FnMut(f64, f64, f64) -> errors::Result<(f64, f64, f64)>,
    {
        let (mut xx, mut yy, mut zz, mut len, stride) = (self.0, self.1, self.2, self.3, self.4);

        if zz.is_null() {
            loop {
                unsafe {
                    f(*xx, *yy, 0.).map_or_else(
                        |_| {
                            *xx = f64::NAN;
                            *yy = f64::NAN;
                        },
                        |(x, y, _)| {
                            *xx = x;
                            *yy = y;
                        },
                    );
                    len -= 1;
                    if len <= 0 {
                        break;
                    }
                    xx = xx.byte_offset(stride);
                    yy = yy.byte_offset(stride);
                }
            }
        } else {
            loop {
                unsafe {
                    f(*xx, *yy, *zz).map_or_else(
                        |_| {
                            *xx = f64::NAN;
                            *yy = f64::NAN;
                            *zz = f64::NAN;
                        },
                        |(x, y, z)| {
                            *xx = x;
                            *yy = y;
                            *zz = z;
                        },
                    );
                    len -= 1;
                    if len <= 0 {
                        break;
                    }
                    xx = xx.byte_offset(stride);
                    yy = yy.byte_offset(stride);
                    zz = zz.byte_offset(stride);
                }
            }
        }
        Ok(())
    }
}

unsafe fn to_degrees(mut x: *mut f64, mut y: *mut f64, mut len: isize, stride: isize) {
    loop {
        *x = (*x).to_degrees();
        *y = (*y).to_degrees();
        len -= 1;
        if len <= 0 {
            break;
        }
        x = x.byte_offset(stride);
        y = y.byte_offset(stride);
    }
}

unsafe fn to_radians(mut x: *mut f64, mut y: *mut f64, mut len: isize, stride: isize) {
    loop {
        *x = (*x).to_radians();
        *y = (*y).to_radians();
        len -= 1;
        if len <= 0 {
            break;
        }
        x = x.byte_offset(stride);
        y = y.byte_offset(stride);
    }
}
