# Changelog

<!--
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/), and this project adheres to [Semantic Versioning](https://semver.org/).
-->

## Unreleased

### Changed


* Use pure rust geodesic implementation for building WASM aeqd in JS context
    - Deprecate 'aeqd' features 
* Use 6th-order series in the third flattening for meridional distance
    - Port of PROJ 9.2 `mlfn` (https://github.com/OSGeo/PROJ/pull/3516)
    - Full double precision; matches PROJ >= 9.2 for `poly`, `sinu`, `eqdc`, `cass`,
      `aeqd` and `tmerc +approx`
* Use series in the third flattening for authalic latitude conversions
    - Port of PROJ 9.7 `latitudes.cpp` (https://github.com/OSGeo/PROJ/pull/4441,
      https://github.com/OSGeo/PROJ/pull/4446)
    - Use `atanh` for computing q (`qsfn`)
    - Matches PROJ >= 9.7 for `aea`, `leac`, `cea`, `laea` and `eqearth`: inverse
      projections are now accurate to double precision (was ~1e-8 deg), forward
      projections differ from PROJ < 9.7 by up to ~5e-8 m


### Added

* Added `cass` Cassini, `robin` Robinson, `eqdc` Equidistant conic 
  and `eqearth` Equal earth projections
    - https://github.com/3liz/proj4rs/pull/59
* Add `tpers/nsper` Tilted/Near-sided perspective projections
    - https://github.com/3liz/proj4rs/pull/58
* Add `sinu` Sinusoidal projection
    - used by MODIS land products and ESRI:54008, 53008 and 102011
    - https://github.com/3liz/proj4rs/pull/57
* Add `omerc` Oblique Mercator
    - https://github.com/3liz/proj4rs/pull/56
* Add `gnom` Gnonomic projection
    - https://github.com/3liz/proj4rs/pull/53
* Add `ortho` Orthographic projection.
    - https://github.com/3liz/proj4rs/pull/52
* Add `poly` American Polyconic projection, used by EPSG:5880, 5530, 29101 and 5472
    - https://github.com/3liz/proj4rs/pull/51
* Add feature for exporting geodesic wasm api

### Fixed

* tmerc spherical: fix inverse latitude sign sign with `lat_0` and forward at the equator
    - https://github.com/3liz/proj4rs/pull/61
* Workaround https://github.com/image-rs/image-tiff/issues/403
    - Enable loading `PROJ-data/ca_nrc/ca_nrc_ntv2_0.tif`
    - Use patched image-tiff branch

## 0.2.1 - 2026-10-09

### Added

* Expose accessors to projection parameters to public api

### Fixed

* Fix standalone nagrids (@)null parameter behavior
* Fix consistent build of nadgrid hierarchy
* Fix trailing (@)null fallback for nadgrids 
* Fix nadgrid subgrid traversal
* Fix UTM zone from `lon_0` on `etmerc`
* Fix `geo_types` transforms
* Fix `krovak` inverse with false easting/northing
* Fix `laea` spherical oblique and south pole forward projection
* Fix forward unit scaling inversion
* Fix spherification formulae
* Support PROJ canonical `k_0`  and aliases
* Fix radians/dms conversion for prime meridians
* Do no skip prime meridian for latlong/geocent

## 0.2.0 - 2026-09-09

### Fixed

* Fix geodetic flattening for ellipsoid built from supplied semi-minor-axis
    - https://github.com/3liz/proj4rs/issues/47
* Fix deprecated field use in `geo_types::Triangle`

### Changed

* Give access to the ellipsoid radius parameters from 'Proj'
* Make module `units` public
* Enforce multi-thread support on native build.

### Added

* Add 'aeqd' projection in default wasm build
* Implement GeoTIFF grids support
* Add `krovak`,`cea` and `mill` in default build.

## 0.1.10 - 2026-03-06

* Bump crs-definitions dependency
    - https://github.com/3liz/proj4rs/pull/42#

* Remove dependency to lazy-static
    - https://github.com/3liz/proj4rs/pull/41

## 0.1.9 - 2025-09-16

### Fixed

* Fix Fix tolerance in inverse grid shift iteration
    - https://github.com/3liz/proj4rs/issues/37

### Changed

* Simplify away some unsafe code
   - https://github.com/3liz/proj4rs/pull/38

## 0.1.8 - 2025-06-18

### Fixed

* Fix NTv2 grid interpolation
    - https://github.com/3liz/proj4rs/issues/32
* Fix Support for wasm-wasip2 target
    - https://github.com/3liz/proj4rs/issues/30
* Fix inverse geos transform and apply more efficient 
  computation: https://github.com/OSGeo/PROJ/pull/4523 

### Added

* DMS notation support for angular proj string values
* Added optional projections: 
    * aeqd 
    * krovak
    * mill
    * cea

## 0.1.7 - 2024-06-10

### Fixed

* Fix axis normalisation

## 0.1.6 - 2024-06-10

### Fixed

* Fix `+nadgrids=@null` as no-op on datum transformation

### Changed 

* Allow 3d inputs in examples/proj4rs

### Added

* Added 'eqc' projection
* Added 'geos' projection 
    - Partially from work from https://github.com/3liz/proj4rs/pull/20

## 0.1.5 - 2024-10-03

### Fixed 

* Fix wrong calculation in laea projection
    - https://github.com/3liz/proj4rs/issues/18

## 0.1.4 - 2024-09-16

### Changed

* Remove wee\_alloc as it's unmaintained
    - https://github.com/3liz/proj4rs/pull/16

## 0.1.3 - 2024-05-18

### Fixed

* fix UB on NodePtr::get
    - https://github.com/3liz/proj4rs/pull/13

### Changed

* Update Vite config to build WASM demo    

## 0.1.2 - 2023-19-11

### Fixed

* Fix `geo-type` feature as optional
    - https://github.com/3liz/proj4rs/pull/11
* Improve documentation
* Fix `Transform` trait signature for WASM
    - https://github.com/3liz/proj4rs/issues/9

### Added

* Add ability to create Projs from EPSG codes
    - https://github.com/3liz/proj4rs/pull/7
* `Transform` implementations.
    - https://github.com/3liz/proj4rs/pull/6
    - Implement for a 2-tuple.
    - Implement for the `geo-types` geometries, them being placed behind a `geo-types` feature flag.

### Changed

* `Transform` trait signature.
    - https://github.com/3liz/proj4rs/pull/6
    - Alias `FnMut(f64, f64, f64) -> Result<(f64, f64, f64)>` behind a `TransformClosure`
    - `transform_coordinates()` takes a mutable reference to `f`, making it easier to layer `Transform` implementations.

## 0.1.1 - 2023-09-07

### Added

* Implement `Clone` on `Proj` type.
    - https://github.com/3liz/proj4rs/pull/2
* Added exemple in README
    - https://github.com/3liz/proj4rs/pull/3
