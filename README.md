[![Crates.io](https://img.shields.io/crates/v/proj4rs)](https://crates.io/crates/proj4rs)
[![Crates.io](https://img.shields.io/crates/d/proj4rs)](https://crates.io/crates/proj4rs)
[![License](https://img.shields.io/badge/license-MIT_OR_Apache_2.0-blue.svg)](https://github.com/image-rs/imageproc/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/Documentation-Published-green)](https://docs.rs/proj4rs/latest/proj4rs/)
[![Demo](https://img.shields.io/badge/Demo-Published-green)](https://docs.3liz.org/proj4rs/)



---

# Rust library for transforming geographic point coordinates from one coordinate system to another.

This is a pure Rust implementation
of the [PROJ.4 project](https://proj.org/en/9.2/faq.html#what-happened-to-proj-4).

The documentation is available on [docs.rs](https://docs.rs/proj4rs/) and the demo on [docs.3liz.org](https://docs.3liz.org/proj4rs/).

The original goal of Proj4rs was to provide a standalone, lightweight, pragmatic implementation of transformations
from CRS to CRS that can be used in Rust and WASM environment.  

C, Python and PHP bindings are provided in order to get reliable and consistent results 
across various contexts (client side JS/server side ....).

## Features and Limitations

Proj4rs is limited to CRS to CRS 2D transformations for an "everyday" use in context where the full
feature set of PROJ is not required or not appropriate. But in all cases, PROJ **is** the source of truth for all 
projections results, and we require proj4rs to stick to this.

Projections are a transcript in Rust of the PROJ projections implementation. The required/expected accuracy
with PROJ 9.8 results is up to 1e-6 m (but in practice it is almost at double precision accuracy).

Have a look at the comparison [projection gallery](https://mdsumner.github.io/proj4rs/).

This port implements the PROJ.4 transformation, which means there's no 3D/4D/orthometric transformation ATM. 
The transformations pipeline use the same WGS84 pivot datum as implemented in PROJ 4.x/5.x.

### Grid supports

Proj4rs support NTV2 and GeoTiff grids.

### Available Projections 

See [PROJECTIONS.md](PROJECTIONS.md) for the list of implemented projections.

## Extra features:

* Geodesic api from `proj4rs-geodesic` crate.

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

Support for WKT is not provided. Instead, there is a dedicated crate for transforming WKT 
strings to proj string.

If you need full support for WKT, please rely on `PROJ` which provides
a great implementation of the standards.

If you want WKT support in WASM, please have a look at:

- https://github.com/3liz/proj4wkt-rs
- https://github.com/frewsxcv/crs-definitions

## Grid shift supports 

Nadgrid is supported for native build and WAMS.

Support for:

* Ntv2 multi grids
* GDAL generated GeoTIFF multi grids

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

You can contribute to this library by going on the [proj4rs](https://github.com/3liz/proj4rs) repository
