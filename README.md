# Proj4rs

This a Proj4 port in Rust. 

This port implemente the Proj4 Api - i.e no 3D/4D/orthometric transformation ATM.

The aim of Proj4rs is to provide at short term the same functionality as the original
proj4 library.

The long term project is to integrate feature from the proj library in its latest
version.

The goal of proj4rs is not to be a remplacement of proj, but instead beeing a light
weight implementation of transformations from crs to crs that could be used
in WASM environment

There is no actual support for WKT, if such supports would exist one day, it would be under
a dedicated crate for transforming proj string to to WKT and vice-versa.

If you need full support for WKT, please rely on proj which provide
a great implementation of the standards.

It is targeted to be WASM compatible for the `wasm32-unknown-unknown` target.

## Grid shift supports 

Currently, only Ntv2 multi grids is supported for native build and WASM.

## Js Api

When compiled for WASM, the library expose a javascript api very similar to proj4js.

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



