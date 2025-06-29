//!
//! Nad grid table
//!
use crate::errors::{Error, Result};
use crate::math::{adjlon, consts::PI};
use crate::transform::Direction;
use std::fmt::{self, Display};

pub(crate) const REL_TOLERANCE_HGRIDSHIFT: f64 = 1.0e-5;

/// Lambda phi pair
#[derive(Debug)]
pub(crate) struct Lp {
    pub lam: f64,
    pub phi: f64,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub(crate) struct GridId([u8; 8]);

impl Default for GridId {
    fn default() -> Self {
        Self::root()
    }
}

impl PartialEq<[u8; 8]> for GridId {
    fn eq(&self, other: &[u8; 8]) -> bool {
        self.0 == *other
    }
}

impl GridId {
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or("<n/a>")
    }

    pub const fn root() -> Self {
        Self([0u8; 8])
    }
}

impl From<[u8; 8]> for GridId {
    fn from(v: [u8; 8]) -> Self {
        Self(v)
    }
}

impl From<u64> for GridId {
    fn from(v: u64) -> Self {
        Self(v.to_ne_bytes())
    }
}

impl From<(u32, u32)> for GridId {
    fn from(p: (u32, u32)) -> Self {
        let mut v: [u8; 8] = [0; 8];
        v[..4].copy_from_slice(&p.0.to_ne_bytes());
        v[4..].copy_from_slice(&p.1.to_ne_bytes());
        Self(v)
    }
}

/// Grid table
#[derive(Debug)]
pub struct Grid {
    pub(crate) id: GridId,
    pub(crate) lineage: GridId,
    pub(crate) ll: Lp,
    pub(crate) ur: Lp,
    pub(crate) del: Lp,
    /// Conversion matrix size
    pub(crate) lim: Lp,
    /// Computed epsilon value
    /// as (fabs(del.lam)+fabs(del.phi))/10000.0
    pub(crate) epsilon: f64,
    /// Conversion matrix: usually stored as f32, f32
    /// and converted to f64, f64
    pub(crate) cvs: Box<[Lp]>,
}

impl Display for Grid {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Id:      {}", self.id.as_str())?;
        writeln!(f, "Lineage: {}", self.lineage.as_str())?;
        writeln!(
            f,
            "Lower Left  ({:>12.7}, {:>12.7})",
            self.ll.lam.to_degrees(),
            self.ll.phi.to_degrees()
        )?;
        writeln!(
            f,
            "Upper Right ({:>12.7}, {:>12.7})",
            self.ur.lam.to_degrees(),
            self.ur.phi.to_degrees()
        )?;
        writeln!(
            f,
            "Size ({} x {})",
            self.lim.lam as usize, self.lim.phi as usize
        )?;
        writeln!(
            f,
            "Delta ({},{})",
            self.del.lam.to_degrees(),
            self.del.phi.to_degrees()
        )
    }
}

impl Grid {
    /// Check if grid is direct child of other.
    #[inline]
    pub fn is_child_of(&self, other: &Grid) -> bool {
        self.lineage == other.id
    }

    #[inline]
    pub fn is_root(&self) -> bool {
        const ROOT: GridId = GridId::root();
        self.lineage == ROOT
    }

    /// Check if the grid match with our point.
    pub fn matches(&self, lam: f64, phi: f64, _z: f64) -> bool {
        !(self.ll.phi - self.epsilon > phi
            || self.ll.lam - self.epsilon > lam
            || self.ll.phi + (self.lim.phi - 1.) * self.del.phi + self.epsilon < phi
            || self.ll.lam + (self.lim.lam - 1.) * self.del.lam + self.epsilon < lam)
    }

    #[inline]
    pub fn num_rows(&self) -> usize {
        self.lim.phi as usize
    }

    #[inline]
    pub fn row_len(&self) -> usize {
        self.lim.lam as usize
    }

    pub fn gs_count(&self) -> usize {
        self.num_rows() * self.row_len()
    }

    pub(crate) fn nad_cvt(
        &self,
        dir: Direction,
        lam: f64,
        phi: f64,
        z: f64,
    ) -> Result<(f64, f64, f64)> {
        match dir {
            Direction::Forward => self.nad_cvt_forward(lam, phi, z),
            Direction::Inverse => self.nad_cvt_inverse(lam, phi, z),
        }
    }

    /// Assume that coordinates matches the grid
    fn nad_cvt_forward(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        let (t_lam, t_phi) = self.nad_intr(
            // normalize input to ll origin
            adjlon(lam - self.ll.lam - PI) + PI,
            phi - self.ll.phi,
        )?;

        Ok((lam + t_lam, phi + t_phi, z))
    }

    fn nad_cvt_inverse(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        const MAX_ITER: usize = 10;
        const TOL: f64 = 1.0e-24;
        const TOL2: f64 = TOL * TOL;

        // normalize input to ll origin
        let (tb_lam, tb_phi) = (adjlon(lam - self.ll.lam - PI) + PI, phi - self.ll.phi);
        let (mut t_lam, mut t_phi) = self.nad_intr(tb_lam, tb_phi)?;

        t_lam = tb_lam - t_lam;
        t_phi = tb_phi - t_phi;

        let mut i = MAX_ITER;
        while i > 0 {
            if let Ok((del_lam, del_phi)) = self.nad_intr(t_lam, t_phi) {
                let (dif_lam, dif_phi) = (t_lam + del_lam - tb_lam, t_phi + del_phi - tb_phi);

                t_lam -= dif_lam;
                t_phi -= dif_phi;

                if dif_lam * dif_lam + dif_phi * dif_phi <= TOL2 {
                    break;
                }

                i -= 1;
            } else {
                // Follows proj5 behavior: returns
                // the first order approximation
                // in case of failure.
                break;
            }
        }

        if i == 0 {
            return Err(Error::InverseGridShiftConvError);
        }

        Ok((adjlon(t_lam + self.ll.lam), t_phi + self.ll.phi, z))
    }

    fn nad_intr(&self, lam: f64, phi: f64) -> Result<(f64, f64)> {
        // Apply bilinear interpolation
        let (t_lam, t_phi) = (lam / self.del.lam, phi / self.del.phi);

        fn _check_lim(t: f64, lim: f64) -> Result<(f64, f64)> {
            let mut i = if t.is_nan() { 0. } else { t.floor() };
            let mut f = t - i;
            if i < 0. {
                if i == -1. && f > 1. - 10. * REL_TOLERANCE_HGRIDSHIFT {
                    i += 1.;
                    f = 0.;
                } else {
                    return Err(Error::PointOutsideNadShiftArea);
                }
            } else {
                let n = i + 1.;
                if n >= lim {
                    if n == lim && f < 10. * REL_TOLERANCE_HGRIDSHIFT {
                        i -= 1.;
                        f = 1.;
                    } else {
                        return Err(Error::PointOutsideNadShiftArea);
                    }
                }
            }
            Ok((i, f))
        }

        let (i_lam, f_lam) = _check_lim(t_lam, self.lim.lam)?;
        let (i_phi, f_phi) = _check_lim(t_phi, self.lim.phi)?;

        let mut index = (i_phi * self.lim.lam + i_lam) as usize;
        let f00 = &self.cvs[index];
        let f10 = &self.cvs[index + 1];
        index += self.lim.lam as usize;
        let f01 = &self.cvs[index];
        let f11 = &self.cvs[index + 1];

        let m00 = (1. - f_lam) * (1. - f_phi);
        let m01 = (1. - f_lam) * f_phi;
        let m10 = f_lam * (1. - f_phi);
        let m11 = f_lam * f_phi;

        Ok((
            m00 * f00.lam + m10 * f10.lam + m01 * f01.lam + m11 * f11.lam,
            m00 * f00.phi + m10 * f10.phi + m01 * f01.phi + m11 * f11.phi,
        ))
    }
}
