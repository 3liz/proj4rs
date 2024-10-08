[![Crates.io](https://img.shields.io/crates/d/proj4rs)](https://crates.io/crates/proj4rs)
[![Documentation](https://img.shields.io/badge/Documentation-Published-green)](https://docs.rs/proj4rs/latest/proj4rs/)
[![Demo](https://img.shields.io/badge/Demo-Published-green)](https://docs.3liz.org/proj4rs/)

---

Rust library for transforming geographic point coordinates
from one coordinate system to another.
This is a pure Rust implementation
of the [PROJ.4 project](https://proj.org/en/9.2/faq.html#what-happened-to-proj-4).

The documentation is available on [docs.rs](https://docs.rs/proj4rs/) and the demo on [docs.3liz.org](https://docs.3liz.org/proj4rs/).

# Features and Limitations

- The aim of Proj4rs is to provide the same functionality as the
[proj4js library](https://github.com/proj4js/proj4js).
- This port implements the PROJ.4 API,
which means there's no 3D/4D/orthometric transformation ATM.
- The goal of Proj4rs is not to be a replacement of PROJ,
but instead being a lightweight implementation of transformations
from CRS to CRS that could be used in Rust and WASM environments.
- This crate does not provide support for WKT. Instead,
there is a dedicated crate for transforming WKT strings to proj string.
- It aims to be WASM compatible for the `wasm32-unknown-unknown` target.
- No installation of external C libraries such as `libproj` or `sqlite3` is needed.

## Basic usage in Rust

Define the coordinate system with proj strings and use the `transform` function.
You can easily get the projection string of any coordinate system
from [EPSG.io](https://epsg.io/).

**Note**: Proj4rs use *radians* as natural angular unit (as does the original proj library)

Example:

```rust
use proj4rs;
use proj4rs::proj::Proj;

// EPSG:5174 - Example
let from = Proj::from_proj_string(concat!(
    "+proj=tmerc +lat_0=38 +lon_0=127.002890277778",
    " +k=1 +x_0=200000 +y_0=500000 +ellps=bessel",
    " +towgs84=-145.907,505.034,685.756,-1.162,2.347,1.592,6.342",
    " +units=m +no_defs +type=crs"
))
.unwrap();

// EPSG:4326 - WGS84, known to us as basic longitude and latitude.
let to = Proj::from_proj_string(concat!(
    "+proj=longlat +ellps=WGS84",
    " +datum=WGS84 +no_defs"
))
.unwrap();

let mut point_3d = (198236.3200000003, 453407.8560000006, 0.0);
proj4rs::transform::transform(&from, &to, &mut point_3d).unwrap();

// Note that WGS84 output from this library is in radians, not degrees.
point_3d.0 = point_3d.0.to_degrees();
point_3d.1 = point_3d.1.to_degrees();

// Output in longitude, latitude, and height.
println!("{} {}",point_3d.0, point_3d.1); // 126.98069676435814, 37.58308534678718
```

## WKT support

If you need full support for WKT, please rely on `proj` which provides
a great implementation of the standards.

If you want WKT support in WASM, please have a look at:

- https://github.com/3liz/proj4wkt-rs
- https://github.com/frewsxcv/crs-definitions

## Grid shift supports 

Nadgrid support is still experimental.
Currently, only Ntv2 multi grids are supported for native build and WASM.

## JavaScript API

When compiled for WASM, the library exposes JavaScript API
that is very similar to that of proj4js.
A thin JavaScript layer provides full compatibility with proj4js
and thus can be used as a proj4js replacement.

Example:

```javascript
let from = new Proj.Projection("+proj=latlong +ellps=GRS80");
let to = new Proj.Projection("+proj=etmerc +ellps=GRS80");
let point = new Proj.Point(2.0, 1.0, 0.0);

// Point is transformed in place
Proj.transform(from, to, point);
```

## Contributing

You can contribute to this library by going on the [proj4rs](./CONTRIBUTING.md) repository
