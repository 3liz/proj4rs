//!
//! Nad grid table
//!
use crate::errors::{Error, Result};
use crate::math::{adjlon, consts::PI};
use crate::transform::Direction;

/// Handle NAD grid and subgrids
#[derive(Debug)]
pub(crate) struct Nadgrid {
    name: String,
    grid: Grid,
    subgrids: Vec<Grid>,
}

impl PartialEq for Nadgrid {
    fn eq(&self, other: &Self) -> bool {
        return self.name == other.name;
    }
}

impl Nadgrid {
    pub(crate) fn new(name: String, grid: Grid, subgrids: Vec<Grid>) -> Self {
        Self {
            name,
            grid,
            subgrids,
        }
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    /// Find the correct grid for an input
    pub(crate) fn find_grid(&self, lam: f64, phi: f64, z: f64) -> Option<&Grid> {
        if self.grid.matches(lam, phi, z) {
            Some(&self.grid)
        } else {
            // Look in sub grids
            self.subgrids.iter().find(|g| g.matches(lam, phi, z))
        }
    }
}

/// Lambda phi pair
#[derive(Debug)]
pub(crate) struct Lp {
    pub(crate) lam: f64,
    pub(crate) phi: f64,
}

/// Grid table
#[derive(Debug)]
pub(crate) struct Grid {
    pub(crate) ll: Lp,
    pub(crate) del: Lp,
    /// Conversion matrix size
    pub(crate) lim: Lp,
    /// Conversion matrix: usually stored as f32, f32
    /// and converted to f64, f64
    pub(crate) flp: Lp,
    /// Computed epsilon value
    /// as (fabs(del.0)+fabs(del.1))/10000.0
    pub(crate) epsilon: f64,
    pub(crate) cvs: Box<[Lp]>,
}

impl Grid {
    /// Check if the grid match with our point.
    pub(crate) fn matches(&self, lam: f64, phi: f64, z: f64) -> bool {
        !(self.ll.phi - self.epsilon > phi
            || self.ll.lam - self.epsilon > lam
            || self.ll.phi + (self.lim.phi - 1.) * self.del.phi + self.epsilon < phi
            || self.ll.lam + (self.lim.lam - 1.) * self.del.lam + self.epsilon < lam)
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

        Ok((lam - t_lam, phi + t_phi, z))
    }

    fn nad_cvt_inverse(&self, lam: f64, phi: f64, z: f64) -> Result<(f64, f64, f64)> {
        const MAX_ITER: usize = 10;
        const TOL: f64 = 1.0e-24; // XXX Check this

        // normalize input to ll origin
        let (tb_lam, tb_phi) = (adjlon(lam - self.ll.lam - PI) + PI, phi - self.ll.phi);
        let (mut t_lam, mut t_phi) = self.nad_intr(tb_lam, tb_phi)?;

        t_lam = tb_lam + t_lam;
        t_phi = tb_phi - t_phi;

        let mut i = MAX_ITER;
        while i > 0 {
            if let Ok((del_lam, del_phi)) = self.nad_intr(t_lam, t_phi) {
                let (diff_lam, diff_phi) = (t_lam - del_lam - tb_lam, t_phi + del_phi - tb_phi);

                if diff_lam * diff_lam + diff_phi * diff_phi <= TOL {
                    break;
                }

                i -= 1;
            } else {
                // Follows proj5 behavior: returns
                // the first order approximation
                // in case of failure.
                i = 0;
                break;
            }
        }

        if i > 0 {
            return Err(Error::InverseGridShiftConvError);
        }

        Ok((adjlon(t_lam + self.ll.lam), t_phi + self.ll.phi, z))
    }

    fn nad_intr(&self, lam: f64, phi: f64) -> Result<(f64, f64)> {
        let (t_lam, t_phi) = (lam / self.del.lam, phi / self.del.phi);

        fn _check_lim(t: f64, lim: f64) -> Result<(f64, f64)> {
            let mut i = t.floor();
            let mut f = t - i;
            if i < 0. {
                if i == -1. && f > 0.99999999999 {
                    i += 1.;
                    f = 0.
                } else {
                    return Err(Error::PointOutsideNadShiftArea);
                }
            } else {
                match i + 1. {
                    n if n == lim && f < 1.0e-11 => {
                        i -= 1.;
                        f = 1.;
                    }
                    n if n > lim => return Err(Error::PointOutsideNadShiftArea),
                    _ => (),
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
