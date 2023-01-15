//!
//! Projection parameters
//!
//!
use crate::errors::{Error, Result};
use std::fmt::Debug;
use std::str::FromStr;

/// Struct holding a pair key/value
pub struct Parameter<'a> {
    pub name: &'a str,
    pub value: Option<&'a str>,
}

// Param -> f64
impl<'a> TryInto<f64> for &Parameter<'a> {
    type Error = Error;

    fn try_into(self) -> Result<f64> {
        self.try_convert::<f64>()
    }
}

// Param -> i32
impl<'a> TryInto<i32> for &Parameter<'a> {
    type Error = Error;

    fn try_into(self) -> Result<i32> {
        self.try_convert::<i32>()
    }
}

// Param -> &str
impl<'a> TryInto<&'a str> for &Parameter<'a> {
    type Error = Error;

    fn try_into(self) -> Result<&'a str> {
        self.value
            .ok_or_else(|| Error::NoValueParameter(self.name.into()))
    }
}

impl<'a> Parameter<'a> {
    pub fn try_convert<F: FromStr>(&self) -> Result<F>
    where
        <F as FromStr>::Err: Debug,
    {
        match self.value.map(F::from_str) {
            None => Err(Error::NoValueParameter(self.name.into())),
            Some(result) => result.map_err(|err| Error::ParameterValueError {
                name: self.name.into(),
                reason: format!("{:?}", err),
            }),
        }
    }
}

/// List of parameters
pub struct ParamList<'a>(Vec<Parameter<'a>>);

impl<'a> ParamList<'a> {
    /// Create a Parameter list from a vector of Params
    pub fn new(params: Vec<Parameter<'a>>) -> Self {
        Self(params)
    }

    /// Return Some(param) if the parameter `name` exists `None` otherwise.
    pub fn get(&self, name: &str) -> Option<&Parameter<'a>> {
        self.0.iter().find(|p| p.name == name)
    }
}
