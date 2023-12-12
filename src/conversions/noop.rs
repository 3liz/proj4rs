//! Reference <https://proj.org/en/9.3/operations/conversions/noop.html>
use crate::*;

#[derive(Debug)]
pub struct NoopConversion;

impl Convert for NoopConversion {
    const NAME: &'static str = "noop";

    type Parameters = ();

    fn new(_: Self::Parameters) -> ProjResult<Self> {
        Ok(Self)
    }

    fn convert(&self, x: f64, y: f64, z: f64) -> ProjResult<(f64, f64, f64)> {
        Ok((x, y, z))
    }
}
