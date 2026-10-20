# proj4rs PHP extension

Native PHP extension for [proj4rs](https://github.com/3liz/proj4rs/), built with
[ext-php-rs](https://github.com/davidcole1340/ext-php-rs).

It exposes the same API as the FFI bindings in `proj4rs-clib/php`
(`Proj4rs\Proj`, `Proj4rs\Transform`, `Proj4rs\Proj4rsException`), without
requiring the `ffi` extension or preloading.

## Building

Requirements: PHP development files (`php-config`) and `libclang`.

```sh
cargo build --release -p proj4rs-php   # => target/release/libproj4rs_php.so
```

The crate is not a workspace default member, so a plain `cargo build` at the
workspace root does not build it.

Add `--features crs-definitions` to enable `EPSG:XXXX` definitions.

The extension is only compatible with the PHP version and thread safety mode
(NTS/ZTS) of the `php-config` it was built with. Set the `PHP_CONFIG`
environment variable to build against another PHP installation.

## Installation

Copy the library into the PHP extension directory (`php-config --extension-dir`)
and enable it in `php.ini`:

```ini
extension=proj4rs.so
```

or try it without installing:

```sh
php -d extension=$PWD/target/release/libproj4rs_php.so script.php
```

## Usage

```php
use Proj4rs\Proj;
use Proj4rs\Transform;

$trans = new Transform('WGS84', '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80');

[$x, $y] = $trans->transform(15.4213696, 47.0766716);
[$xs, $ys] = $trans->transform([15.4213696, 2.35], [47.0766716, 48.85]);
$points = $trans->transformPoints([[15.4213696, 47.0766716], [2.35, 48.85]]);
```

See `proj4rs-clib/php/README.md` for the full API.

## Tasks

```sh
cargo make php.test     # Run the PHP test suite against the extension
cargo make php.package  # Build target/php/proj4rs-php-ext-<version>-php<X.Y>-<nts|zts>-<target>.tar.gz
```
