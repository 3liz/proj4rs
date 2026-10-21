# Projection gallery

A visual and numerical check of every proj4rs projection against PROJ,
run in a browser on the wasm build.

Each panel is drawn entirely by proj4rs: the Natural Earth 110m coastline
(grey) and a graticule (light blue) are projected point by point through
the proj4js-compatible wrapper (`js/proj4.js`). The red rings are the
graticule intersections projected by PROJ. Where the two agree the rings
sit on the line crossings. Under each map the page reports, over all
intersections:

- forward: largest distance between proj4rs and PROJ, in projected units
- inverse: largest error of the proj4rs inverse of the PROJ coordinates, in degrees
- visibility mismatches: points one library projects and the other refuses
- inverse failures: PROJ coordinates the proj4rs inverse cannot invert

Verdicts: `identical` (bit for bit), `matches` (forward within `tol`,
default 1e-3, and inverse within 1e-7 deg), `known` (a difference
recorded in `projections.json` with a note), `differs` (a new difference).

## Files

- `projections.json`: the catalogue. One entry per map: `id`, `proj`
  string, and optionally `lon`/`lat` window (default whole world), `step`
  for graticule intersections (default 30 deg), `tol`, `note`, and
  `expect: "differs"` for a known difference.
- `make_reference.py`: computes the PROJ intersections (through pyproj,
  `Transformer.from_pipeline` on the bare proj string, so no CRS
  interpretation and no datum step) into `reference.json`, and downloads
  the coastline into `data/`.
- `index.html`, `gallery.js`: the page.
- `render.py`: opens the page in headless Chromium, writes
  `out/gallery.png`, one PNG per projection and `out/results.json`, prints
  a table, and exits 1 on any unexpected result.

The lon/lat source system for each map takes the earth model of the
projection (`+ellps`, `+R`, `+a`, `+b`, ...), so no datum shift is involved.

## Running

```bash
# 1. wasm module into js/pkg
cd proj4rs && cargo make wasm && cd ..

# 2. PROJ reference and coastline (needs pyproj)
python3 js/gallery/make_reference.py

# 3a. look at it
cd js && python3 -m http.server   # then open http://localhost:8000/gallery/

# 3b. or render and check headlessly (needs the playwright python package)
python3 js/gallery/render.py      # --only ortho,gnom for a subset
```

## Adding a projection

Add an entry to `projections.json`, rerun `make_reference.py`, then
`render.py`. Pick a window where the projection is defined (transverse
and conic projections need a lon/lat window, azimuthal ones usually the
whole world). For a documented difference add `note` and
`expect: "differs"` so the check stays green until it is fixed.
