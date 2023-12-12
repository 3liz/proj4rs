//! Refernce: <https://proj.org/en/9.3/operations/conversions/axisswap.html>

use crate::*;

#[derive(Debug, Clone)]
pub struct AxisswapConversion {
    ordering: AxisswapOrdering,
}

impl Convert for AxisswapConversion {
    const NAME: &'static str = "axisswap";
    type Parameters = AxisswapOrdering;

    fn new(ordering: Self::Parameters) -> ProjResult<Self> {
        Ok(Self { ordering })
    }

    fn convert(&self, x: f64, y: f64, z: f64) -> ProjResult<(f64, f64, f64)> {
        let output = self.ordering.apply_ordering([x, y, z]);
        Ok((output[0], output[1], output[2]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disallows_missing_order_parameter() {
        assert!(Conversion::from_proj_string("+proj=axisswap").is_err(),)
    }

    #[test]
    fn converts_from_proj_string() {
        let conversion = Conversion::from_proj_string("+proj=axisswap +order=2,1").unwrap();
        let mut points = (1., 2., 0.);
        conversion.convert(&mut points).unwrap();
        assert_eq!((2., 1., 0.), points);
    }
}

pub use ordering::AxisswapOrdering;
mod ordering {
    use std::str::FromStr;

    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct Flip(bool);

    #[derive(Debug, Clone)]
    #[cfg_attr(test, derive(PartialEq))]
    pub struct AxisswapOrdering([(u8, Flip); 3]);

    impl ConvertParameters for AxisswapOrdering {
        fn from_parameter_list(parameter_list: &ParamList) -> ProjResult<Self> {
            parameter_list
                .get("order")
                .ok_or(ProjError::NoValueParameter)?
                .value
                .ok_or(ProjError::NoValueParameter)?
                .parse::<AxisswapOrdering>()
        }
    }

    impl AxisswapOrdering {
        const AXIS_COUNT: usize = 3;

        pub fn apply_ordering(&self, input: [f64; Self::AXIS_COUNT]) -> [f64; Self::AXIS_COUNT] {
            let mut output = [0.; 3];

            (0..Self::AXIS_COUNT).for_each(|input_index| {
                let (final_location, flip) = self
                    .0
                    .iter()
                    .enumerate()
                    .find_map(|(final_location, (axis_number, flip))| {
                        (axis_number == &(input_index as u8)).then_some((final_location, flip))
                    })
                    .expect("no axis number in order array");

                output[final_location] = {
                    let mut value = input[input_index];

                    if flip == &Flip(true) {
                        value *= -1.0;
                    }

                    value
                };
            });

            output
        }
    }

    impl FromStr for AxisswapOrdering {
        type Err = ProjError;

        fn from_str(ordering_str: &str) -> Result<Self, Self::Err> {
            let mut found_axes: [Option<(u8, Flip)>; 3] = [None; Self::AXIS_COUNT];

            for (found_axis_index, value_str) in ordering_str.split(',').enumerate() {
                let value = value_str.parse::<i8>().map_err(|_| {
                    ProjError::InvalidParameterValue(
                        "unable to parse comma separated value into integer",
                    )
                })?;

                let (axis_number, flip) = match value < 0 {
                    true => (-value as u8, Flip(true)),
                    false => (value as u8, Flip(false)),
                };

                if axis_number == 0 {
                    return Err(ProjError::InvalidParameterValue(
                        "axis value out of range: 0",
                    ));
                }

                if axis_number > Self::AXIS_COUNT as u8 {
                    return Err(ProjError::InvalidParameterValue(
                        "axis value larger than number of dimensions",
                    ));
                }

                found_axes[found_axis_index] = Some((axis_number, flip));
            }

            let mut to_swap: [Option<(u8, Flip)>; 3] = [None; Self::AXIS_COUNT];

            // fill unspecifed values in to_swap with no_op
            for maybe_unspecified_axis in 1..(Self::AXIS_COUNT + 1) {
                if !found_axes.iter().any(|element| {
                    element
                        .is_some_and(|(found_axis, _)| found_axis == maybe_unspecified_axis as u8)
                }) {
                    let unspecified_axis = maybe_unspecified_axis;
                    to_swap[unspecified_axis - 1] = Some((unspecified_axis as u8, Flip(false)));
                }
            }

            // fill found axis in whichever locations are not already occupied
            for (found_axis, flip) in found_axes.into_iter().flatten() {
                let Some(unoccupied_value) = to_swap.iter_mut().find(|element| element.is_none())
                else {
                    return Err(ProjError::InvalidParameterValue(
                        "duplicate axes are disallowed",
                    ));
                };

                *unoccupied_value = Some((found_axis, flip));
            }

            let mut ordering = [(0, Flip(false)); Self::AXIS_COUNT];

            to_swap
                .into_iter()
                // all positions should now be filled
                .map(|element| element.expect("to swap location not specified"))
                .enumerate()
                // decrement by one to represent index locations
                .for_each(|(index, (to_swap_axis, flip))| {
                    ordering[index] = (to_swap_axis - 1, flip)
                });

            Ok(Self(ordering))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        impl AxisswapOrdering {
            fn mock(order: [u8; Self::AXIS_COUNT]) -> Self {
                Self([
                    (order[0], Flip(false)),
                    (order[1], Flip(false)),
                    (order[2], Flip(false)),
                ])
            }
        }

        #[test]
        fn performs_order_swap() {
            // All possible permutations (3! = 6)
            assert_eq!(
                AxisswapOrdering::mock([0, 1, 2]).apply_ordering([1., 2., 3.]),
                [1., 2., 3.]
            );

            assert_eq!(
                AxisswapOrdering::mock([0, 2, 1]).apply_ordering([1., 2., 3.]),
                [1., 3., 2.]
            );

            assert_eq!(
                AxisswapOrdering::mock([1, 0, 2]).apply_ordering([1., 2., 3.]),
                [2., 1., 3.]
            );

            assert_eq!(
                AxisswapOrdering::mock([1, 2, 0]).apply_ordering([1., 2., 3.]),
                [2., 3., 1.]
            );

            assert_eq!(
                AxisswapOrdering::mock([2, 0, 1]).apply_ordering([1., 2., 3.]),
                [3., 1., 2.]
            );

            assert_eq!(
                AxisswapOrdering::mock([2, 1, 0]).apply_ordering([1., 2., 3.]),
                [3., 2., 1.]
            );
        }

        #[test]
        fn performs_axis_flip() {
            assert_eq!(
                AxisswapOrdering([(0, Flip(true)), (1, Flip(true)), (2, Flip(false))])
                    .apply_ordering([1., -2., 3.]),
                [-1., 2., 3.]
            );
        }

        #[test]
        fn parses_valid_order() {
            assert_eq!(
                AxisswapOrdering::mock([2, 0, 1]),
                "3,1,2".parse::<AxisswapOrdering>().unwrap()
            )
        }

        #[test]
        fn parses_only_necessary_pair() {
            assert_eq!(
                AxisswapOrdering::mock([1, 0, 2]),
                "2,1".parse::<AxisswapOrdering>().unwrap()
            );
            assert_eq!(
                AxisswapOrdering::mock([2, 1, 0]),
                "3,1".parse::<AxisswapOrdering>().unwrap()
            );
            assert_eq!(
                AxisswapOrdering::mock([0, 2, 1]),
                "3,2".parse::<AxisswapOrdering>().unwrap()
            );
        }

        #[test]
        fn parses_singular_value() {
            assert_eq!(
                AxisswapOrdering::mock([0, 1, 2]),
                "3".parse::<AxisswapOrdering>().unwrap()
            )
        }

        #[test]
        fn parses_direction() {
            assert_eq!(
                AxisswapOrdering([(0, Flip(true)), (1, Flip(false)), (2, Flip(true))]),
                "-1,2,-3".parse::<AxisswapOrdering>().unwrap()
            )
        }

        #[test]
        fn disallows_axis_zero() {
            assert!("1,0,2".parse::<AxisswapOrdering>().is_err())
        }

        #[test]
        fn disallows_axis_greater_than_total_count() {
            assert!("1,5,3".parse::<AxisswapOrdering>().is_err())
        }

        #[test]
        fn disallows_duplicates() {
            assert!("1,2,1".parse::<AxisswapOrdering>().is_err())
        }
    }
}
