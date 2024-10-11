
The documentation is available on [docs.rs](https://docs.rs/proj4rs/) and the demo on [docs.3liz.org](https://docs.3liz.org/proj4rs/).

# Build locally a package

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

This will create a web package in js/pkg

## Build for npm

```bash
cargo make wasm_bundle
```

This will create a npm bundler package in js/pkg-bundler


Packages are created in the *js/*  directory at the root of the repository.
