# Projections implemented in proj4rs

Projections are selected with the `+proj=<name>` parameter. Aliases share the
implementation of their parent projection.

| Name | Description | PROJ documentation |
|------|-------------|--------------------|
| `latlong`, `longlat` | Geographic coordinates (no projection) | [latlon](https://proj.org/en/stable/operations/conversions/latlon.html), [lonlat](https://proj.org/en/stable/operations/conversions/lonlat.html) |
| `geocent`, `cart` | Geocentric (cartesian) coordinates | [cart](https://proj.org/en/stable/operations/conversions/cart.html) |
| `aea` | Albers Equal Area (conic, sph & ell) | [aea](https://proj.org/en/stable/operations/projections/aea.html) |
| `leac` | Lambert Equal Area Conic (alias of `aea`) | [leac](https://proj.org/en/stable/operations/projections/leac.html) |
| `aeqd` | Azimuthal Equidistant (azimuthal, sph & ell) | [aeqd](https://proj.org/en/stable/operations/projections/aeqd.html) |
| `cass` | Cassini (cylindrical, sph & ell) | [cass](https://proj.org/en/stable/operations/projections/cass.html) |
| `cea` | Equal Area Cylindrical (cylindrical, sph & ell) | [cea](https://proj.org/en/stable/operations/projections/cea.html) |
| `eqc` | Equidistant Cylindrical (Plate Carrée) | [eqc](https://proj.org/en/stable/operations/projections/eqc.html) |
| `eqdc` | Equidistant Conic (conic, sph & ell) | [eqdc](https://proj.org/en/stable/operations/projections/eqdc.html) |
| `eqearth` | Equal Earth (pseudocylindrical, sph & ell) | [eqearth](https://proj.org/en/stable/operations/projections/eqearth.html) |
| `etmerc` | Extended Transverse Mercator (Poder/Engsager "exact" algorithm) | [tmerc](https://proj.org/en/stable/operations/projections/tmerc.html) |
| `utm` | Universal Transverse Mercator (based on `etmerc`) | [utm](https://proj.org/en/stable/operations/projections/utm.html) |
| `geos` | Geostationary Satellite View | [geos](https://proj.org/en/stable/operations/projections/geos.html) |
| `gnom` | Gnomonic (azimuthal, sph & ell) | [gnom](https://proj.org/en/stable/operations/projections/gnom.html) |
| `krovak` | Krovak (oblique conformal conic, Czech/Slovak systems) | [krovak](https://proj.org/en/stable/operations/projections/krovak.html) |
| `labrd` | Laborde (cylindrical, ellipsoidal, for Madagascar) | [labrd](https://proj.org/en/stable/operations/projections/labrd.html) |
| `laea` | Lambert Azimuthal Equal Area (azimuthal, sph & ell) | [laea](https://proj.org/en/stable/operations/projections/laea.html) |
| `lcc` | Lambert Conformal Conic (conic, sph & ell) | [lcc](https://proj.org/en/stable/operations/projections/lcc.html) |
| `merc` | Mercator (cylindrical, sph & ell) | [merc](https://proj.org/en/stable/operations/projections/merc.html) |
| `webmerc` | Web Mercator / Pseudo Mercator (alias of `merc`) | [webmerc](https://proj.org/en/stable/operations/projections/webmerc.html) |
| `mill` | Miller Cylindrical (spherical) | [mill](https://proj.org/en/stable/operations/projections/mill.html) |
| `moll` | Mollweide (pseudocylindrical, spherical) | [moll](https://proj.org/en/stable/operations/projections/moll.html) |
| `wag4` | Wagner IV (pseudocylindrical, variant of `moll`) | [wag4](https://proj.org/en/stable/operations/projections/wag4.html) |
| `wag5` | Wagner V (pseudocylindrical, variant of `moll`) | [wag5](https://proj.org/en/stable/operations/projections/wag5.html) |
| `omerc` | Oblique Mercator (cylindrical, sph & ell) | [omerc](https://proj.org/en/stable/operations/projections/omerc.html) |
| `ortho` | Orthographic (azimuthal, sph & ell) | [ortho](https://proj.org/en/stable/operations/projections/ortho.html) |
| `poly` | American Polyconic (sph & ell) | [poly](https://proj.org/en/stable/operations/projections/poly.html) |
| `robin` | Robinson (pseudocylindrical, spherical) | [robin](https://proj.org/en/stable/operations/projections/robin.html) |
| `sinu` | Sinusoidal / Sanson-Flamsteed (pseudocylindrical, sph & ell) | [sinu](https://proj.org/en/stable/operations/projections/sinu.html) |
| `somerc` | Swiss Oblique Mercator (cylindrical, ellipsoidal, for CH1903) | [somerc](https://proj.org/en/stable/operations/projections/somerc.html) |
| `stere` | Stereographic (azimuthal, sph & ell) | [stere](https://proj.org/en/stable/operations/projections/stere.html) |
| `ups` | Universal Polar Stereographic (alias of `stere`) | [ups](https://proj.org/en/stable/operations/projections/ups.html) |
| `sterea` | Oblique Stereographic Alternative (azimuthal, sph & ell) | [sterea](https://proj.org/en/stable/operations/projections/sterea.html) |
| `tmerc` | Transverse Mercator (Poder/Engsager by default, Evenden/Snyder for the spherical case or with `+approx`) | [tmerc](https://proj.org/en/stable/operations/projections/tmerc.html) |
| `tpers` | Tilted Perspective (azimuthal, spherical) | [tpers](https://proj.org/en/stable/operations/projections/tpers.html) |
| `nsper` | Near-sided Perspective (alias of `tpers`) | [nsper](https://proj.org/en/stable/operations/projections/nsper.html) |
