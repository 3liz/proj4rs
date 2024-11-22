//!
//! Crate errors
//!

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("{0}")]
    InputStringError(&'static str),
    #[error("Missing value for parameter {0}")]
    NoValueParameter(String),
    #[error("Cannot retrieve value for parameter")]
    ParameterValueError(String),
    #[error("Missing projection name")]
    MissingProjectionError,
    #[error("Unrecognized datum")]
    InvalidDatum,
    #[error("Unrecognized ellipsoid")]
    InvalidEllipsoid,
    #[error("{0}")]
    InvalidParameterValue(&'static str),
    #[error("Invalid coordinate dimension")]
    InvalidCoordinateDimension,
    #[error("Latitude out of range")]
    LatitudeOutOfRange,
    #[error("NAD grid not available")]
    NadGridNotAvailable,
    #[error("Parent grid not found")]
    NadGridParentNotFound,
    #[error("Inverse grid shift failed to converge.")]
    InverseGridShiftConvError,
    #[error("Point outside of NAD outside Shift area")]
    PointOutsideNadShiftArea,
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
    #[error("Invalid number of coordinates")]
    InvalidNumberOfCoordinates,
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
    #[error("Invalid Ntv2 grid format: {0}")]
    InvalidNtv2GridFormat(&'static str),
    #[error("IO error")]
    IoError(#[from] std::io::Error),
    #[error("UTF8 error")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("Grid file not found {0}")]
    GridFileNotFound(String),
    #[error("Unknown grid format")]
    UnknownGridFormat,
    #[error("Numerical argument too  large")]
    ArgumentTooLarge,
}

pub type Result<T, E = Error> = std::result::Result<T, E>;
