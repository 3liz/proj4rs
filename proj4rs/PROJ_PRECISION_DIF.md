# PROJ 9.8.1 precision against exact values

Date: 2026-11-16. PROJ 9.8.1 through pyproj 3.8.0 (local `.venv`), proj4rs `main` with the
changes described in `PRECISIONS.md`.

## Method

* **Points:** the cases of `proj4rs/tests/precision_tests.rs`. They were chosen near the poles,
  the equator and the projection centres, where PROJ loses precision, so they are not
  representative of typical accuracy.
* **Exact values:** 50- to 60-digit computations of the closed formulas. For the ellipsoidal
  round-trip cases, the exact value is the original (lon, lat).
* **Units:** forward errors are distances in metres between the projected coordinates; inverse
  errors are great-circle distances in metres (R = 6371000 m for the sphere, a = 6378137 m for
  GRS80).

## PROJ 9.8.1 errors, worst first

| Rank | Case | PROJ 9.8.1 error |
|---|---|---|
| 1 | `laea` sphere, `+lat_0=45`, inverse at (90, -10): returns longitude 0 instead of 90 | **9,815 km** |
| 2 | `eqearth` sphere, inverse at the pole: `MAX_Y` 1.3e-14 too small | 0.74 m |
| 3 | `moll` sphere, forward at (-179, -89.99999): Newton iteration gives up | 0.145 m |
| 4 | `laea` GRS80 polar, forward at 1e-6° from the pole: returns (0, 0) | 0.112 m |
| 5 | `laea` / `stere` GRS80 oblique (`lat_0=-30`), round trip at 1e-7° from the pole | 8.3 cm |
| 6 | `stere` GRS80 oblique (`lat_0=45`), round trip at 1e-6° from the pole | 2.2 cm |
| 7 | `laea` GRS80 polar and oblique, inverse near the poles | 1.1–1.7 cm |
| 8 | `gnom` sphere, inverse close to 90° from the centre | 0.3–1.6 cm |
| 9 | `stere`, `ortho`, `aeqd`, `laea`, `tpers` spherical inverses, and `tmerc` sphere at 1e-7° from the equator | ~1.1 cm |
| 10 | `moll` inverse at (-179, -89.99999) | 4.8 mm |
| 11 | `laea` polar forward at 1e-5° from the pole | 0.45 mm |
| 12 | `laea` and `aeqd` near the projection centre | 0.22–0.25 mm |
| 13 | `moll` at (100, 89.99): forward / inverse | 13 µm / 2 µm |

Only 3 of the 38 points have PROJ within 1e-6 m of the exact value in both directions: the polar
`gnom` point and the `gnom` and `tpers` centres.

## Comments

* **proj4rs at the same points:** the worst error is 3.4e-7 m, in the `moll` inverse at
  (-179, -89.99999). Apart from the two `moll` points, every point is within 1.9e-9 m, and most
  errors are about 1.4e-9 m, a single rounding step in the latitude near the poles
  (2.2e-16 rad). That is about 10⁵ times smaller than PROJ's centimetre-level errors.
* **Remaining `moll` error in proj4rs:** near the pole, Mollweide amplifies tiny absolute errors
  (x grows like the colatitude to the power 2/3).
  * Forward: `FRAC_PI_2` is 6.1e-17 below π/2, which gives about 1.4e-7 m in x at 1e-5° from
    the pole. Adding its low part (`6.123233995736766e-17`) to the colatitude would remove it.
  * Inverse: rounding y divided by the radius and storing √2 as a double, inside the
    difference `c_y - |y|` (about 7.7e-10 there), give the 3.3e-7 m. Correcting √2 with its
    low part would halve that; the first rounding happens in the transform pipeline before the
    projection is called.
* **Test reference:** the `moll` (-179, -89.99999) reference in `precision_tests.rs` was computed
  from the exact decimal degrees, not from the double radians proj4rs receives. For that input
  the exact value is x = -591.0783918757571, not -591.078391808469. The 6.7e-8 m difference
  doesn't change any result in this table.
* **Ordinary points:** away from these places, proj4rs and PROJ 9.8.1 agree within 1e-6 m at
  all 4,600 points of `tests/proj_parity_tests.rs`, on a sphere and on GRS80. Those points have
  no exact values, so they show agreement, not accuracy.
