// Statements to load wasm module in nodejs REPL
let Proj
import("../../pkg-node/proj4rs.js").then(module => { Proj = module });
