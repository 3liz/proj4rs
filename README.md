# Proj4rs

This a Proj4 port in Rust. 

This port implemente the Proj4 Api - i.e no 3D/4D/orthometric transformation ATM.

The aim of Proj4rs is to provide the same functionality as the original
proj4js library.

The goal of proj4rs is not to be a remplacement of proj, but instead beeing a light
weight implementation of transformations from crs to crs that could be used
in WASM environment

This crate does not provide support for WKT, instead, there is a dedicated crate for transforming 
WKT strings to proj string.

It is targeted to be WASM compatible for the `wasm32-unknown-unknown` target.

Documentation on [doc.rs](https://docs.rs/proj4rs/)

## WKT support

If you need full support for WKT, please rely on `proj` which provides
a great implementation of the standards.

If you want WKT support in WASM please have a look at https://github.com/3liz/proj4wkt-rs

## Grid shift supports 

Currently, only Ntv2 multi grids is supported for native build and WASM.

## Js Api

When compiled for WASM, the library expose a javascript api very similar to proj4js. A thin
javascript layer provide full compatibility with proj4js and thus can be used as a proj4js
replacement.

Example:

```javascript
let from = new Proj.Projection("+proj=latlong +ellps=GRS80");
let to = new Proj.Projection("+proj=etmerc +ellps=GRS80");
let point = new Proj.Point(2.0, 1.0, 0.0);

// Point is transformed in place
Proj.transform(from, to, point);
```

## Compiling for WASM

Install [wasm-pack](https://rustwasm.github.io/wasm-pack/book/)

```
wasm-pack build --target web --no-default-features
```

Or if you have installed [cargo-make](https://sagiegurari.github.io/cargo-make/), use the following
command:

```
cargo make wasm
```

### Running the WASM example

There is a [`index.html`] file for testing the WASM module in a navigator.

For security reason you need to run it from a server; you can pop up 
a server from python with the following command:

```
python3 -m http.server
```
 
The server will automatically serve the `index.html` file in the current directory.



