use crate::*;

#[derive(Debug)]
pub enum Conversion {
    Axisswap(AxisswapConversion),
    Noop(NoopConversion),
}

impl Conversion {
    pub fn from_proj_string(proj_str: &str) -> ProjResult<Self> {
        let parameter_list = projstring::parse(proj_str)?;

        let conversion_name = parameter_list
            .get("proj")
            .and_then(|parameter| parameter.value)
            .ok_or(ProjError::MissingProjectionError)?;

        match conversion_name {
            AxisswapConversion::NAME => {
                AxisswapConversion::from_params_list(&parameter_list).map(Self::Axisswap)
            }
            NoopConversion::NAME => {
                NoopConversion::from_params_list(&parameter_list).map(Self::Noop)
            }
            _ => Err(ProjError::InvalidParameterValue("unrecognized projection")),
        }
    }

    pub fn convert<T: Transform>(&self, points: &mut T) -> ProjResult<()> {
        match self {
            Conversion::Axisswap(conversion) => {
                points.transform_coordinates(&mut |x, y, z| conversion.convert(x, y, z))
            }
            Conversion::Noop(conversion) => {
                points.transform_coordinates(&mut |x, y, z| conversion.convert(x, y, z))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converts_from_proj_str() {
        let mut points = (1.0, 2.0);
        Conversion::from_proj_string("+proj=noop")
            .unwrap()
            .convert(&mut points)
            .unwrap();
    }
}
