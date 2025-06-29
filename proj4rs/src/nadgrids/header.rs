//!
//! Nadgrid parser
//!
use crate::errors::{Error, Result};
use crate::nadgrids::grid::GridId;
use std::io::Read;

#[derive(Copy, Clone)]
pub(crate) enum Endianness {
    Be = 0,
    Le = 1,
}

#[cfg(target_endian = "big")]
impl Endianness {
    pub fn native() -> Self {
        Endianness::Be
    }
    pub fn other() -> Self {
        Endianness::Le
    }
}
#[cfg(target_endian = "little")]
impl Endianness {
    pub fn native() -> Self {
        Endianness::Le
    }
    pub fn other() -> Self {
        Endianness::Be
    }
}

/// Generic header struct
pub(crate) struct Header<const N: usize> {
    buf: [u8; N],
    pub endian: Endianness,
}

impl<const N: usize> Header<N> {
    pub fn new() -> Self {
        Self::new_endian(Endianness::native())
    }

    pub fn new_endian(endian: Endianness) -> Self {
        Self {
            buf: [0u8; N],
            endian,
        }
    }

    pub fn rebind<const M: usize>(&self) -> Header<M> {
        Header::<M>::new_endian(self.endian)
    }

    #[inline]
    pub fn read<R: Read>(&mut self, read: &mut R) -> Result<&Self> {
        read.read_exact(&mut self.buf)?;
        Ok(self)
    }

    #[inline]
    pub fn read_partial<R: Read>(&mut self, read: &mut R) -> Result<usize> {
        read.read(&mut self.buf).map_err(Error::from)
    }

    pub fn get_f64(&self, offset: usize) -> f64 {
        match self.endian {
            Endianness::Be => f64::from_be_bytes(self.buf[offset..offset + 8].try_into().unwrap()),
            Endianness::Le => f64::from_le_bytes(self.buf[offset..offset + 8].try_into().unwrap()),
        }
    }

    pub fn get_f32(&self, offset: usize) -> f32 {
        match self.endian {
            Endianness::Be => f32::from_be_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
            Endianness::Le => f32::from_le_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
        }
    }

    pub fn get_u32(&self, offset: usize) -> u32 {
        match self.endian {
            Endianness::Be => u32::from_be_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
            Endianness::Le => u32::from_le_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
        }
    }

    #[inline]
    pub fn get_str(&self, offset: usize, len: usize) -> Result<&str> {
        std::str::from_utf8(&self.buf[offset..offset + len]).map_err(Error::from)
    }

    #[inline]
    pub fn get_u8(&self, offset: usize) -> u8 {
        self.buf[offset]
    }

    pub fn get_id(&self, offset: usize) -> GridId {
        let mut buf = [0u8; 8];
        buf.copy_from_slice(&self.buf[offset..offset + 8]);
        buf.into()
    }

    pub fn cmp_str(&self, offset: usize, s: &str) -> bool {
        self.get_str(offset, s.len())
            .map(|x| x == s)
            .unwrap_or(false)
    }
}

pub mod error_str {
    pub const ERR_INVALID_HEADER: &str = "Invalid header";
    pub const ERR_GSCOUNT_NOT_MATCHING: &str = "GS COUNT not matching";
}
