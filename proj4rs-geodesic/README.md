Geodesic routines
=================

This is a pure Rust port of the C implementation of the geodesic
algorithms of [GeographicLib](https://geographiclib.sourceforge.io/) as
implemented in the [Proj library](https://proj.org/en/stable/).

Works on wasm32.

It provides:

* the solution of the direct and inverse geodesic problems (`Geodesic`),
* the computation of points along a geodesic (`GeodesicLine`),
* the perimeter and area of geodesic polygons (`Polygon`).

Results are identical to those of the reference C implementation.

This crate is used as an optional dependency of the [proj4rs](https://lib.rs/crates/proj4rs)
and was created initially for the support of the [Azimuthal Equidistant (aeqd)](https://proj.org/en/stable/operations/projections/aeqd.html) projection.


## Motivations

The C math library does not link *as-is* in JS WASM context, and thus, projections based 
on this lib were missing in the WASM/JS context.  Rewriting the library in pure idiomatic
rust solve the problem.
