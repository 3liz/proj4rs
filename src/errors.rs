//!
//! Crate errors
//!

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("{0}")]
    InputStringError(&'static str),
    #[error("No value for parameter")]
    NoValueParameter,
    #[error("Cannot retrieve value for parameter")]
    ParameterValueError,
    #[error("Missing projection name")]
    MissingProjectionError,
    #[error("Unrecognized datum")]
    InvalidDatum,
    #[error("Unrecognized ellipsoid")]
    InvalidEllipsoid,
    #[error("{0}")]
    InvalidParameterValue(&'static str),
    #[error("Latitude out of range")]
    LatitudeOutOfRange,
    #[error("NAD grid not available")]
    NoNADGridAvailable,
    #[error("Invalid 'towgs84' string")]
    InvalidToWGS84String,
    #[error("Invalid axis")]
    InvalidAxis,
    #[error("Unrecognized format")]
    UnrecognizedFormat,
    #[error("Latitude or longitude over range")]
    LatOrLongExceedLimit,
    #[error("Nan value for coordinate")]
    NanCoordinateValue,
    #[error("Coordinate out of range")]
    CoordinateOutOfRange,
    #[error("Projection not found")]
    ProjectionNotFound,
    #[error("No forward projection defined for dest projection")]
    NoForwardProjectionDefined,
    #[error("No inverse projection defined for src projection")]
    NoInverseProjectionDefined,
    #[error("ProjErrConicLatEqual")]
    ProjErrConicLatEqual,
    #[error("Tolerance condition not satisfied")]
    ToleranceConditionError,
    #[error("Non convergence of phi2 calculation")]
    NonInvPhi2Convergence,
    #[error("Failed no compute forward projection")]
    ForwardProjectionFailure,
    #[error("Failed no compute inverse projection")]
    InverseProjectionFailure,
    #[error("Invalid UTM zone")]
    InvalidUtmZone,
    #[error("An ellipsoid is required")]
    EllipsoidRequired,
    #[error("Coordinate transform outside projection domain")]
    CoordTransOutsideProjectionDomain,
    #[error("No convergence for inv. meridian distance")]
    InvMeridDistConvError,
    #[error("JS parse error")]
    JsParseError,
}

pub type Result<T> = std::result::Result<T, Error>;
