
# Notes about precision

Proj4rs use PROJ transformations (`proj`, `cs2cs`) results as refererence for all tests using PROJ 9.4.0 to
PROJ latest. Precision required is usually 10^-6 (up to micrometers) for length/angles.

Is some case the precision may be less (or more) which depends of the PROJ version outputs you are comparing to, but never
more than 10^-4.

## Differences with PROJ

### 1. Meridian distance

- **proj4rs** (src/math/mlfn.rs) still uses the formula PROJ had before 9.2. It is Snyder's equation 3-21 as 
  a series in e² up to e⁸. The inverse (`inv_mlfn`) runs a Newton iteration on that series until the step 
  is below 1e-11 rad.
- **PROJ 9.2.0** replaced it (PR #3516 (https://github.com/OSGeo/PROJ/pull/3516), Karney). 
  The new formula is a 6th-order series in the third flattening n, and the inverse is a direct series with 
  no iteration. Series in n converge much faster than series in e², so PROJ now reaches full double precision 
  for |f| ≤ 1/150.
- **PROJ 9.7.0** moved mlfn.cpp onto the new auxiliary-latitude code (#4446 (https://github.com/OSGeo/PROJ/pull/4446)). 
  This changed nothing visible at 1e-9 m.


### 2. Authalic latitude (eqearth, and also aea, cea, laea)

- **proj4rs and PROJ ≤ 9.6:** the forward direction uses qsfn, a closed form with a log term. 
  The inverse (authset/authlat) is a 3-term series in e² from Adams (1921), truncated at e⁶.
- **PROJ 9.7.0 changed both directions:**
  - #4441 (https://github.com/OSGeo/PROJ/pull/4441) added an "exact" inverse solved by Newton's method.
    The forward q now uses atanh instead of log, which keeps relative accuracy near the equator. 
    This came from an issue about ETRS89-LAEA accuracy.
  - #4446 (https://github.com/OSGeo/PROJ/pull/4446) switched to series in n for small flattening,
    evaluated with Clenshaw summation.

The effect on eqearth WGS84 with the 9.8.1 test values:

- **Forward**: about 1e-8 m different. This is only last-digit noise from atanh versus log and from the reformulation.
- **Inverse**: the old Adams series is off by up to ~1.2e-8 deg (~1.3 mm).
  For example, the point at −30° comes back as −29.999999988 in both proj4rs and PROJ 9.4.0, while 9.8.1 returns −30. 
  That's why the eqearth inverse test needs 1e-7 deg.

In proj4rs vs PROJ < 9.7.0 Authalic forward/inverse are identical

## Math libraries

Some Rust math functions are implemented in Rust's std itself rather than calling platform library: 
`asinh`, `acosh`, `atanh` and `ln_1p`  (on some targets) and their results can differ highly from `glibc`'s.  

Some of theses or redefined in PROJ (`ln_1p`), while others are redefined in Proj4rs 
(asinh in `proj4rs/src/math.rs`, `remquo`, `remainder`, `atanh` `in proj4rs-geodesic/src/math.rs`) - see
the comments in their respective files for more details.

The difference between platform math libraries are around 1e-9 m, which is far beyond the requested precision.




