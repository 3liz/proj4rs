//! Wasm

/// Use [cargo-run-wasm](https://github.com/rukai/cargo-run-wasm) to build an example for web
///
/// Usage:
/// ```
/// cargo run-wasm --package [example_name]
/// ```
/// Generally:
/// ```
/// cargo run-wasm -p proj4rs
/// ```
fn main() {
    cargo_run_wasm::run_wasm_cli_with_css("body { margin: 0px; }");
}
