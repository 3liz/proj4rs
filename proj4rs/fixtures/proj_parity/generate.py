"""Generate proj_9_8_1.csv: PROJ outputs at ordinary points, for tests/proj_parity_tests.rs.

Run with pyproj built on PROJ 9.8.1 (e.g. the repository .venv):

    .venv/bin/python -I proj4rs/fixtures/proj_parity/generate.py

Each row is `def;lon;lat;x;y;inv_lon;inv_lat`: (x, y) is PROJ's forward of (lon, lat), in
degrees, and (inv_lon, inv_lat) is PROJ's inverse of that (x, y).

The points avoid the places where PROJ itself loses precision (see commit 45ff983,  73340c4, 2fadbc3 and
tests/precision_tests.rs): latitudes within 1 degree of the poles, the antipode and horizon of the
azimuthal projections, the map outline of moll, and, for tmerc, longitudes from 85 to 95 degrees
from the central meridian (PROJ 9.8.1 does not include OSGeo/PROJ#4945, whose easting loses digits
close to +/-90 degrees) and latitudes within 1 degree of the equator beyond 95 degrees (where its
northing loses digits).

Regional projections (somerc, labrd) are sampled in a longitude/latitude box around their area of
use, on the figure of their definition.
"""

import math
import pathlib
import random
import sys

import pyproj
from pyproj import Transformer

if pyproj.proj_version_str != "9.8.1":
    sys.exit(f"PROJ 9.8.1 is required, found {pyproj.proj_version_str}")

POINTS = 100
SPHERE = "+R=6371000"
GRS80 = "+ellps=GRS80"

# (definition without the figure, figures, max angular distance from (lon_0, lat_0) in degrees,
#  max |lon| in degrees[, min |lon| in degrees, min |lat| in degrees])
SETUPS = [
    ("+proj=laea +lat_0=90", (SPHERE, GRS80), 170, 180),
    ("+proj=laea +lat_0=-90", (SPHERE, GRS80), 170, 180),
    ("+proj=laea +lat_0=45", (SPHERE, GRS80), 170, 180),
    ("+proj=laea +lat_0=0", (SPHERE, GRS80), 170, 180),
    ("+proj=laea +lat_0=-30", (SPHERE, GRS80), 170, 180),
    ("+proj=stere +lat_0=90 +k_0=0.994", (SPHERE, GRS80), 150, 180),
    ("+proj=stere +lat_0=45 +k_0=1", (SPHERE, GRS80), 150, 180),
    ("+proj=stere +lat_0=0 +k_0=1", (SPHERE, GRS80), 150, 180),
    ("+proj=stere +lat_0=-30 +k_0=0.9996", (SPHERE, GRS80), 150, 180),
    ("+proj=ortho +lat_0=90", (SPHERE, GRS80), 85, 180),
    ("+proj=ortho +lat_0=45", (SPHERE, GRS80), 85, 180),
    ("+proj=ortho +lat_0=0", (SPHERE, GRS80), 85, 180),
    ("+proj=aeqd +lat_0=90", (SPHERE, GRS80), 170, 180),
    ("+proj=aeqd +lat_0=45", (SPHERE, GRS80), 170, 180),
    ("+proj=aeqd +lat_0=0", (SPHERE, GRS80), 170, 180),
    ("+proj=aeqd +lat_0=-30", (SPHERE, GRS80), 170, 180),
    ("+proj=gnom +lat_0=90", (SPHERE,), 80, 180),
    ("+proj=gnom +lat_0=45", (SPHERE,), 80, 180),
    ("+proj=gnom +lat_0=0", (SPHERE,), 80, 180),
    ("+proj=tpers +lat_0=45 +h=3000000", (SPHERE,), 42, 180),
    ("+proj=moll", (SPHERE,), 360, 179.9),
    ("+proj=eqearth", (SPHERE, GRS80), 360, 180),
    ("+proj=tmerc +k=0.999", (SPHERE,), 360, 85),
    # spherical tmerc beyond 90 degrees from the central meridian
    ("+proj=tmerc +k=0.999", (SPHERE,), 360, 180, 95, 1),
    ("+proj=tmerc +k=0.9996 +lat_0=40", (SPHERE,), 360, 180, 95, 1),
]

# Regional projections: (full definition, (min lon, max lon), (min lat, max lat)) in degrees.
# The geographic CRS uses the figure (+R= or +ellps=) of the definition.
REGIONS = [
    # Swiss LV95 (EPSG:2056), Bessel
    (
        "+proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1"
        " +x_0=2600000 +y_0=1200000 +ellps=bessel",
        (5.5, 10.8),
        (45.6, 48.0),
    ),
    ("+proj=somerc +lat_0=-30 +lon_0=20 +k_0=0.9996 +x_0=0 +y_0=0 +ellps=GRS80", (10, 30), (-40, -20)),
    ("+proj=somerc +lat_0=46 +lon_0=8 +x_0=0 +y_0=0 +R=6371000", (0, 16), (40, 52)),
    # Madagascar Laborde (EPSG:29702), International 1924
    (
        "+proj=labrd +lat_0=-18.9 +lon_0=46.4372291666667 +azi=18.9 +k_0=0.9995"
        " +x_0=400000 +y_0=800000 +ellps=intl",
        (43, 51),
        (-26, -12),
    ),
]


def figure(definition):
    return " ".join(p for p in definition.split() if p.startswith(("+R=", "+ellps=")))


def lat_0(definition):
    for p in definition.split():
        if p.startswith("+lat_0="):
            return float(p.split("=")[1])
    return 0.0


def angular_distance(lon, lat, lat0):
    p1, p2, dl = math.radians(lat0), math.radians(lat), math.radians(lon)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def main():
    random.seed(20261116)
    out = pathlib.Path(__file__).with_name("proj_9_8_1.csv")
    rows = []
    for base, figures, max_dist, max_lon, *mins in SETUPS:
        min_lon, min_lat = mins or (0.0, 0.0)
        for fig in figures:
            definition = f"{base} {fig} +x_0=0 +y_0=0 +lon_0=0"
            t = Transformer.from_crs(
                f"+proj=longlat {fig}", definition + " +type=crs", always_xy=True
            )
            n = 0
            while n < POINTS:
                lon = random.uniform(-max_lon, max_lon)
                lat = random.uniform(-89.0, 89.0)
                if abs(lon) < min_lon or abs(lat) < min_lat:
                    continue
                if angular_distance(lon, lat, lat_0(base)) > max_dist:
                    continue
                x, y = t.transform(lon, lat)
                if not (math.isfinite(x) and math.isfinite(y)):
                    continue
                inv_lon, inv_lat = t.transform(x, y, direction="INVERSE")
                if not (math.isfinite(inv_lon) and math.isfinite(inv_lat)):
                    continue
                rows.append(f"{definition};{lon!r};{lat!r};{x!r};{y!r};{inv_lon!r};{inv_lat!r}")
                n += 1
    for definition, (lon_min, lon_max), (lat_min, lat_max) in REGIONS:
        t = Transformer.from_crs(
            f"+proj=longlat {figure(definition)}", definition + " +type=crs", always_xy=True
        )
        for _ in range(POINTS):
            lon = random.uniform(lon_min, lon_max)
            lat = random.uniform(lat_min, lat_max)
            x, y = t.transform(lon, lat)
            inv_lon, inv_lat = t.transform(x, y, direction="INVERSE")
            rows.append(f"{definition};{lon!r};{lat!r};{x!r};{y!r};{inv_lon!r};{inv_lat!r}")
    out.write_text(
        f"# PROJ {pyproj.proj_version_str} outputs, generated by generate.py\n"
        "# def;lon;lat;x;y;inv_lon;inv_lat\n" + "\n".join(rows) + "\n"
    )
    print(f"{len(rows)} rows written to {out}")


if __name__ == "__main__":
    main()
