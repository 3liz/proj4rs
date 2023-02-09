//!
//! Nadgrid parser
//!
use crate::errors::{Error, Result};
use crate::nadgrids::GridId;
use std::io::Read;

#[derive(Copy, Clone)]
enum Endianness {
    Be = 0,
    Le = 1,
}

#[cfg(target_endian = "big")]
impl Endianness {
    fn native() -> Self {
        Endianness::Be
    }
    fn other() -> Self {
        Endianness::Le
    }
}
#[cfg(target_endian = "little")]
impl Endianness {
    fn native() -> Self {
        Endianness::Le
    }
    fn other() -> Self {
        Endianness::Be
    }
}

/// Generic header struct
struct Header<const N: usize> {
    buf: [u8; N],
    endian: Endianness,
}

impl<const N: usize> Header<N> {
    fn new() -> Self {
        Self::new_endian(Endianness::native())
    }

    fn new_endian(endian: Endianness) -> Self {
        Self {
            buf: [0u8; N],
            endian,
        }
    }

    fn rebind<const M: usize>(&self) -> Header<M> {
        Header::<M>::new_endian(self.endian)
    }

    #[inline]
    fn read<R: Read>(&mut self, read: &mut R) -> Result<&Self> {
        read.read_exact(&mut self.buf)?;
        Ok(self)
    }

    fn get_f64(&self, offset: usize) -> f64 {
        match self.endian {
            Endianness::Be => f64::from_be_bytes(self.buf[offset..offset + 8].try_into().unwrap()),
            Endianness::Le => f64::from_le_bytes(self.buf[offset..offset + 8].try_into().unwrap()),
        }
    }

    fn get_f32(&self, offset: usize) -> f32 {
        match self.endian {
            Endianness::Be => f32::from_be_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
            Endianness::Le => f32::from_le_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
        }
    }

    fn get_u32(&self, offset: usize) -> u32 {
        match self.endian {
            Endianness::Be => u32::from_be_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
            Endianness::Le => u32::from_le_bytes(self.buf[offset..offset + 4].try_into().unwrap()),
        }
    }

    #[inline]
    fn get_str(&self, offset: usize, len: usize) -> Result<&str> {
        std::str::from_utf8(&self.buf[offset..offset + len]).map_err(Error::from)
    }

    #[inline]
    fn get_u8(&self, offset: usize) -> u8 {
        self.buf[offset]
    }

    #[inline]
    fn get_id(&self, offset: usize) -> GridId {
        let mut buf = [0u8; 8];
        buf.copy_from_slice(&self.buf[offset..offset + 8]);
        buf.into()
    }
}

mod error_str {
    pub(super) const ERR_INVALID_HEADER: &str = "Invalid header";
    pub(super) const ERR_GSCOUNT_NOT_MATCHING: &str = "GS COUNT not matching";
}

// Parsers
mod ntv2;
