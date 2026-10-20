# proj4rs PHP bindings

PHP bindings for [proj4rs](https://github.com/3liz/proj4rs/) using the
[FFI](https://www.php.net/manual/en/book.ffi.php) extension over the `proj4rs_c`
C library.

## Requirements

* PHP >= 8.1 with the `ffi` extension
* The `libproj4rs_c` shared library:

```sh
cargo build --release -p proj4rs-clib   # => target/release/libproj4rs_c.so
```

Add `--features crs-definitions` to enable `EPSG:XXXX` definitions.

## Loading the library

The library is resolved in this order:

1. An explicit call to `Proj4rs\Library::load('/path/to/libproj4rs_c.so')`
2. The `PROJ4RS` FFI scope, if `proj4rs.h` was preloaded (see below)
3. The `PROJ4RS_LIB` environment variable
4. The library bundled in the package `lib/` directory
5. The platform library name (`libproj4rs_c.so`, `libproj4rs_c.dylib`, `proj4rs_c.dll`),
   resolved by the dynamic linker

With the default `ffi.enable=preload` setting, FFI is only allowed from the CLI
or from preloaded code. For a web server (php-fpm), preload the header from
your `opcache.preload` script:

```php
<?php
FFI::load('/path/to/proj4rs-clib/php/proj4rs.h');
```

Edit `FFI_LIB` at the top of `proj4rs.h` if the library is not in the linker search path.

## Usage

```php
use Proj4rs\Proj;
use Proj4rs\Transform;

$src = new Proj('WGS84');
$dst = new Proj('+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80');

$trans = new Transform($src, $dst);

// Scalars
[$x, $y] = $trans->transform(15.4213696, 47.0766716);   // 4732659.007426, 2677630.726961

// Arrays of coordinates
[$xs, $ys] = $trans->transform([15.4213696, 2.35], [47.0766716, 48.85]);

// List of points
$points = $trans->transformPoints([[15.4213696, 47.0766716], [2.35, 48.85, 0.0]]);
```

* If `convert` is `true` (the default), latlong coordinates are in degrees.
* Points that fail to project are set to `NAN`.
* Errors raise `Proj4rs\Proj4rsException`.

`Proj` exposes `projname()`, `isLatlong()`, `isGeocent()`, `axis()`,
`isNormalizedAxis()`, `toMeter()` and `units()`.

## Tests

```sh
PROJ4RS_LIB=../../target/release/libproj4rs_c.so php tests/test_all.php
```
