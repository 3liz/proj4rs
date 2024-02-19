use crate::*;

pub trait Convert: Sized {
    const NAME: &'static str;

    type Parameters: ConvertParameters;

    fn new(parameters: Self::Parameters) -> ProjResult<Self>;

    fn convert(&self, x: f64, y: f64, z: f64) -> ProjResult<(f64, f64, f64)>;

    fn from_params_list(parameter_list: &ParamList) -> ProjResult<Self> {
        Self::new(<Self::Parameters as ConvertParameters>::from_parameter_list(parameter_list)?)
    }
}

pub trait ConvertParameters: Sized {
    fn from_parameter_list(parameter_list: &ParamList) -> ProjResult<Self>;
}

impl ConvertParameters for () {
    fn from_parameter_list(_: &ParamList) -> ProjResult<Self> {
        Ok(())
    }
}
