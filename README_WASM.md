
The documentation is available on [docs.rs](https://docs.rs/proj4rs/) and the demo on [docs.3liz.org](https://docs.3liz.org/proj4rs/).

## Compiling for WASM

Install [wasm-pack](https://rustwasm.github.io/wasm-pack/book/)

```bash
wasm-pack build --target web --no-default-features
```

Or if you have installed [cargo-make](https://sagiegurari.github.io/cargo-make/), use the following
command:

```bash
cargo make wasm
```

### Running the WASM example

There is a [`index.html`] file for testing the WASM module in a navigator.

For security reasons, you need to run it from a server.
You can start a Python server with the following command:

```bash
python3 -m http.server
```

The server will automatically serve the `index.html` file in the current directory.


## Build for npm

```
cargo make wasm_bundle
```

This will create a npm bundler package in pkg-bundler
