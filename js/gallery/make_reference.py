#!/usr/bin/env python3
"""Build the PROJ reference data for the proj4rs projection gallery.

For every entry in projections.json, project the graticule intersections
(every `step` degrees inside the entry's lon/lat window) with PROJ through
pyproj, and write them to reference.json. The gallery page draws these as
red rings over the proj4rs rendering and reports the difference.

PROJ is called with Transformer.from_pipeline on the bare proj string, so
the projection runs exactly as written (no CRS interpretation, no datum
step). The Natural Earth 110m coastline is downloaded to data/ if missing.

Usage: python3 make_reference.py
Needs: pyproj (any PROJ >= 9), network access for the first run.
"""
import json
import math
import os
import urllib.request

import pyproj
from pyproj import Transformer
from pyproj.enums import TransformDirection

HERE = os.path.dirname(os.path.abspath(__file__))
COAST_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "master/geojson/ne_110m_coastline.geojson"
)
DEFAULTS = {"lon": [-180, 180], "lat": [-90, 90], "step": 30}


def fetch_coastline():
    path = os.path.join(HERE, "data", "ne_110m_coastline.geojson")
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        print("downloading", COAST_URL)
        urllib.request.urlretrieve(COAST_URL, path)


def frange(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [lo + i * step for i in range(n + 1)]


def finite(v):
    return math.isfinite(v) and abs(v) < 1e30


def reference(entry):
    e = {**DEFAULTS, **entry}
    lons = frange(e["lon"][0], e["lon"][1], e["step"])
    lats = frange(e["lat"][0], e["lat"][1], e["step"])
    is_longlat = e["proj"].split()[0] in ("+proj=longlat", "+proj=latlong")
    if not is_longlat:
        t = Transformer.from_pipeline(
            "+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad "
            "+step " + e["proj"]
        )
    points = []
    for lat in lats:
        for lon in lons:
            if is_longlat:
                x, y = lon, lat
            else:
                x, y = t.transform(lon, lat, direction=TransformDirection.FORWARD)
            points.append([lon, lat, x, y] if finite(x) and finite(y) else [lon, lat, None, None])
    return points


def main():
    fetch_coastline()
    entries = json.load(open(os.path.join(HERE, "projections.json")))
    out = {
        "proj_version": pyproj.proj_version_str,
        "pyproj_version": pyproj.__version__,
        "defaults": DEFAULTS,
        "entries": {e["id"]: reference(e) for e in entries},
    }
    with open(os.path.join(HERE, "reference.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))
    n = sum(len(v) for v in out["entries"].values())
    print(f"wrote reference.json: {len(entries)} projections, {n} points, PROJ {pyproj.proj_version_str}")


if __name__ == "__main__":
    main()
