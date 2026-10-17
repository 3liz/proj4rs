//! Low level math helpers.
//!
//! These are transcriptions of the static helpers of `geodesic.c`,
//! together with the C99 functions (`remainder`, `remquo`) that are not
//! available in the Rust standard library.

/// Quarter turn in degrees
pub(crate) const QD: f64 = 90.0;
/// Half turn in degrees
pub(crate) const HD: f64 = 2.0 * QD;
/// Full turn in degrees
pub(crate) const TD: f64 = 2.0 * HD;

/// Radians per degree
pub(crate) const DEGREE: f64 = std::f64::consts::PI / HD;

/// Exact power of 2 for exponents in the normal range
const fn pow2(e: i64) -> f64 {
    f64::from_bits(((1023 + e) as u64) << 52)
}

/// `sqrt(f64::MIN_POSITIVE)`
pub(crate) const TINY: f64 = pow2(-511);
pub(crate) const TOL0: f64 = f64::EPSILON;
/// Increase multiplier in defn of tol1 from 100 to 200 to fix inverse case
/// 52.784459512564 0 -52.784459512563990912 179.634407464943777557
/// which otherwise failed for Visual Studio 10 (Release and Debug)
pub(crate) const TOL1: f64 = 200.0 * TOL0;
/// `sqrt(f64::EPSILON)`
pub(crate) const TOL2: f64 = pow2(-26);
/// Check on bisection interval
pub(crate) const TOLB: f64 = TOL0;
pub(crate) const XTHRESH: f64 = 1000.0 * TOL2;

pub(crate) const MAXIT1: u32 = 20;
pub(crate) const MAXIT2: u32 = MAXIT1 + f64::MANTISSA_DIGITS + 10;

#[inline]
pub(crate) fn sq(x: f64) -> f64 {
    x * x
}

/// Error-free sum: `u + v = s + t = round(u + v) + t`
///
/// Returns `(s, t)`.
#[inline]
pub(crate) fn sumx(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    let up = up - u;
    let vpp = vpp - v;
    let t = if s != 0.0 { 0.0 - (up + vpp) } else { s };
    (s, t)
}

/// Evaluate polynomial `p[0] * x^n + ... + p[n]` with `n = p.len() - 1`
#[inline]
pub(crate) fn polyval(p: &[f64], x: f64) -> f64 {
    match p.split_first() {
        Some((&y, rest)) => rest.iter().fold(y, |y, &c| y * x + c),
        None => 0.0,
    }
}

/// IEEE 754 remainder: `x - n * y` where `n` is `x / y` rounded to the
/// nearest integer (ties to even). The result is exact.
pub(crate) fn remainder(x: f64, y: f64) -> f64 {
    let p = y.abs();
    // Now |r| < 2p (exact)
    let r = x % (2.0 * p);
    let mut a = r.abs();
    let half = 0.5 * p;
    // Sterbenz lemma: the subtractions below are exact
    if a > half {
        a -= p;
        if a >= half {
            a -= p;
        }
    }
    // Zero result takes the sign of x
    if r.is_sign_negative() { -a } else { a }
}

/// C99 `remquo` returning the remainder and the low order bits
/// (at least 2) of the quotient.
pub(crate) fn remquo(x: f64, y: f64) -> (f64, i32) {
    // x = k * 4y + r4, the quotient of x / y modulo 4 is the
    // quotient of r4 / y modulo 4 (ties are resolved the same since
    // 4k is even).
    let r4 = x % (4.0 * y);
    let r = remainder(r4, y);
    // r4 - r is an exact multiple of y
    let q = ((r4 - r) / y).round();
    (r, if q.is_nan() { 0 } else { q as i32 })
}

/// Inverse hyperbolic tangent
///
/// `f64::atanh` from std is less accurate than the C99 `atanh` used by the
/// reference implementation, so use the fdlibm algorithm.
pub(crate) fn atanh(x: f64) -> f64 {
    let y = x.abs();
    let t = if y < 0.5 {
        let t = y + y;
        0.5 * (t + t * y / (1.0 - y)).ln_1p()
    } else {
        0.5 * ((y + y) / (1.0 - y)).ln_1p()
    };
    t.copysign(x)
}

#[inline]
pub(crate) fn norm2(sinx: f64, cosx: f64) -> (f64, f64) {
    let r = sinx.hypot(cosx);
    (sinx / r, cosx / r)
}

pub(crate) fn ang_normalize(x: f64) -> f64 {
    let y = remainder(x, TD);
    if y.abs() == HD { HD.copysign(x) } else { y }
}

#[inline]
pub(crate) fn lat_fix(x: f64) -> f64 {
    if x.abs() > QD { f64::NAN } else { x }
}

/// Compute `y - x` reduced to `[-180, 180]`, accurately.
///
/// Returns `(d, e)` where `d + e` is the exact difference.
pub(crate) fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    // Use remainder instead of AngNormalize, since we treat boundary cases
    // later taking account of the error
    let (d, t) = sumx(remainder(-x, TD), remainder(y, TD));
    // This second sum can only change d if abs(d) < 128, so don't need to
    // apply remainder yet again.
    let (mut d, t) = sumx(remainder(d, TD), t);
    // Fix the sign if d = -180, 0, 180.
    if d == 0.0 || d.abs() == HD {
        // If t == 0, take sign from y - x
        // else (t != 0, implies d = +/-180), d and t must have opposite signs
        d = d.copysign(if t == 0.0 { y - x } else { -t });
    }
    (d, t)
}

/// Round an angle so that small values are coarsened to multiples of 1/16
/// (this reduces the chance of loss of precision with tiny values).
pub(crate) fn ang_round(x: f64) -> f64 {
    const Z: f64 = 1.0 / 16.0;
    let y = x.abs();
    let w = Z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    let y = if w > 0.0 { Z - w } else { y };
    y.copysign(x)
}

#[inline]
fn sincos_quadrant(r: f64, q: i32, x: f64) -> (f64, f64) {
    let (s, c) = r.sin_cos();
    let (sinx, cosx) = match (q as u32) & 3 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    };
    // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
    // special values from F.10.1.12
    let cosx = cosx + 0.0;
    // special values from F.10.1.13
    let sinx = if sinx == 0.0 { sinx.copysign(x) } else { sinx };
    (sinx, cosx)
}

/// Sine and cosine of an angle in degrees.
///
/// In order to minimize round-off errors, this function exactly reduces
/// the argument to the range [-45, 45] before converting it to radians.
pub(crate) fn sincosd(x: f64) -> (f64, f64) {
    let (r, q) = remquo(x, QD);
    // now abs(r) <= 45
    sincos_quadrant(r * DEGREE, q, x)
}

/// Sine and cosine of `x + t` in degrees, where `t` is a small correction.
pub(crate) fn sincosde(x: f64, t: f64) -> (f64, f64) {
    let (r, q) = remquo(x, QD);
    let r = ang_round(r + t);
    // now abs(r) <= 45
    sincos_quadrant(r * DEGREE, q, x)
}

/// `atan2` with the result in degrees.
///
/// In order to minimize round-off errors, this function rearranges the
/// arguments so that result of atan2 is in the range [-pi/4, pi/4] before
/// converting it to degrees and mapping the result to the correct
/// quadrant.
pub(crate) fn atan2d(mut y: f64, mut x: f64) -> f64 {
    let mut q = 0;
    if y.abs() > x.abs() {
        std::mem::swap(&mut x, &mut y);
        q = 2;
    }
    if x.is_sign_negative() {
        x = -x;
        q += 1;
    }
    // here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
    let ang = y.atan2(x) / DEGREE;
    match q {
        1 => HD.copysign(y) - ang,
        2 => QD - ang,
        3 => -QD + ang,
        _ => ang,
    }
}

/// Accumulator with twice the normal precision
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct Accumulator {
    s: f64,
    t: f64,
}

impl Accumulator {
    pub fn sum(&self) -> f64 {
        self.s
    }

    pub fn add(&mut self, y: f64) {
        let (z, u) = sumx(y, self.t);
        let (s, t) = sumx(z, self.s);
        if s == 0.0 {
            self.s = u;
            self.t = t;
        } else {
            self.s = s;
            self.t = t + u;
        }
    }

    /// Return accumulator + y (but don't add to accumulator).
    pub fn sum_with(&self, y: f64) -> f64 {
        let mut a = *self;
        a.add(y);
        a.s
    }

    pub fn negate(&mut self) {
        self.s = -self.s;
        self.t = -self.t;
    }

    /// Reduce to [-y/2, y/2].
    pub fn remainder(&mut self, y: f64) {
        self.s = remainder(self.s, y);
        self.add(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remainder() {
        assert_eq!(remainder(5.0, 2.0), 1.0);
        assert_eq!(remainder(7.0, 2.0), -1.0);
        assert_eq!(remainder(3.0, 2.0), -1.0);
        assert_eq!(remainder(-3.0, 2.0), 1.0);
        assert_eq!(remainder(540.0, TD), -180.0);
        assert_eq!(remainder(-540.0, TD), 180.0);
        assert_eq!(remainder(180.0, TD), 180.0);
        assert_eq!(remainder(190.0, TD), -170.0);
        assert!(remainder(-720.0, TD).is_sign_negative());
        assert!(remainder(720.0, TD).is_sign_positive());
        assert!(remainder(f64::INFINITY, TD).is_nan());
        assert!(remainder(f64::NAN, TD).is_nan());
    }

    // Compare with the C99 implementation from the system libm
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    #[test]
    fn test_remainder_libm() {
        unsafe extern "C" {
            #[link_name = "remainder"]
            fn c_remainder(x: f64, y: f64) -> f64;
            #[link_name = "remquo"]
            fn c_remquo(x: f64, y: f64, q: *mut i32) -> f64;
        }
        // Simple LCG for reproducible pseudo random values
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let area0 = 4.0 * std::f64::consts::PI * 40_589_732_499_314.76;
        for i in 0..200_000 {
            let scale = [1.0, 90.0, 720.0, 1e6, 1e15, 1e300][i % 6];
            let x = (2.0 * next() - 1.0) * scale;
            // Exact multiples of 45 to check ties
            let x = if i % 7 == 0 {
                (x / 45.0).round() * 45.0
            } else {
                x
            };
            for y in [QD, TD, 2.0 * TD, area0] {
                let r = remainder(x, y);
                let c = unsafe { c_remainder(x, y) };
                assert_eq!(r.to_bits(), c.to_bits(), "remainder({x}, {y})");
            }
            let (r, q) = remquo(x, QD);
            let mut cq = 0;
            let c = unsafe { c_remquo(x, QD, &mut cq) };
            assert_eq!(r.to_bits(), c.to_bits(), "remquo({x})");
            assert_eq!(q & 3, cq & 3, "remquo({x}) quotient");
        }
    }

    #[test]
    fn test_remquo() {
        assert_eq!(remquo(0.0, QD), (0.0, 0));
        assert_eq!(remquo(45.0, QD), (45.0, 0));
        assert_eq!(remquo(135.0, QD), (-45.0, 2));
        assert_eq!(remquo(225.0, QD), (45.0, 2));
        assert_eq!(remquo(-100.0, QD), (-10.0, -1));
        assert_eq!(remquo(270.0, QD).0, 0.0);
        assert_eq!(remquo(270.0, QD).1 & 3, 3);
        assert_eq!(remquo(-270.0, QD).1 & 3, 1);
    }

    #[test]
    fn test_accumulator() {
        let mut a = Accumulator::default();
        a.add(1.0);
        a.add(1e-20);
        a.add(-1.0);
        assert_eq!(a.sum(), 1e-20);
    }
}
