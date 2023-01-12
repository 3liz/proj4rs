//!
//! # Proj string parser
//!
//! ## Projection string grammar
//!
//! ```text
//! <projstring> :: =
//!   +<ident>=<value> { <projstring> }... |
//!   +<ident> { <projstring> }...
//!
//! <ident> ::= [_0-9a_Z]+
//!
//! <value> ::= " { <word> }... " | <word>
//!
//! <word> ::= [^\s]+
//! ```
//!
//! Possible parameters of a projection may be
//!
//! ## Cartograpic projection parameters:
//!
//! see <https://proj.org/usage/projections.html#cartographic-projection>
//!
//! +a      : semi major axis
//! +axis   : Axis orientation
//! +b      : semi minor axis
//! +ellps  : ellipsoid name
//! +k      : Scaling factor (deprecated)
//! +k_0    : Scaling factor
//! +lat_0  : Latitude of origin
//! +long_0 : Central meridian
//! +pm     : Alternate prime meridian
//! +proj   : Projection name
//! +units  : meters, US survey feet, etc.
//! +x_0    : False easting
//! +y_0    : False northing
//!
//!
//! ## Geodetic transformations parameters
//!
//! see <https://proj.org/usage/transformation.html#geodetic-transformation>
//! +datum    : Datum name
//! +to_meter : Multiplier to convert map units to 1.0m
//! +towgs84  : 3 or 7 term datum transform parameters
//! +nadgrids : Filename of NTv2 grid file to use for datum transforms
//!
//!
//! ## Ellipsoid parameters
//!
//! see <https://proj.org/usage/ellipsoids.html>
//!
//! +rf    : reverse flattening
//! +a      : semi major axis
//! +b      : semi minor axis
//!
//!
//! ## Per projections parameters  
//!
//! These parameters depends on the projection used.
//! One must refer to the projection definition.
//!
use crate::errors::{Error, Result};
use std::ops::ControlFlow;

struct Parser {}

impl Parser {
    /// Parse parameter name as valid identifier
    ///
    /// Return error if the string is not a valid
    /// identifier: i.e not [0-9a-zA-Z_]+
    fn parse_identifier(s: &str) -> Result<(&str, &str)> {
        // Get the identifiant
        let rv = s.chars().try_fold(Ok(0usize), |len, c| {
            if c.is_whitespace() || c == '=' {
                ControlFlow::Break(len)
            } else if !c.is_alphanumeric() && c != '_' {
                // Invalid character for identifier
                ControlFlow::Break(Err(Error::InputStringError("Invalid parameter name")))
            } else {
                ControlFlow::Continue(len.map(|len| len + c.len_utf8()))
            }
        });

        match rv {
            ControlFlow::Break(res) => res.map(|len| (&s[..len], &s[len..])),
            ControlFlow::Continue(_) => Ok((s, "")),
        }
    }

    /// Get the next quoted or unquoted token from the input string
    fn unquote_next(s: &str) -> Result<(&str, &str)> {
        let s = s.trim_start();
        if s.starts_with('\"') {
            // Check if string part is terminated by a quote,
            // if not, continue with the next part.
            // Inner quotes not separated by a whitespace is left
            // as part of the token.
            let s = s.split_once('\"').unwrap().1;
            match s
                .split_inclusive(|c: char| c.is_whitespace())
                .try_fold(0usize, |len, s| {
                    let offset = s.len();
                    let s = s.trim_end();
                    if s.ends_with('\"') {
                        ControlFlow::Break(len + s.len() - 1)
                    } else {
                        ControlFlow::Continue(len + offset)
                    }
                }) {
                ControlFlow::Break(len) => Ok((&s[..len], &s[(len + 1)..])),
                ControlFlow::Continue(len) => {
                    Err(Error::InputStringError("Unterminated quoted string"))
                }
            }
        } else {
            Ok(s.split_once(|c: char| c.is_whitespace()).unwrap_or((s, "")))
        }
    }

    /// Returns the first token from the input str
    fn token(s: &str) -> Result<(&str, Option<&str>, &str)> {
        let s = s.trim_start();
        if s.is_empty() {
            Ok(("", None, ""))
        } else if s.starts_with('+') {
            let (_, rest) = s.split_once('+').unwrap(); // Swallow '+'

            let (name, rest) = Self::parse_identifier(rest)?;
            if name.is_empty() {
                Err(Error::InputStringError("Empty parameter name"))
            } else {
                let rest = rest.trim_start();
                if rest.starts_with('=') {
                    let (value, rest) = Self::unquote_next(rest.split_once('=').unwrap().1)?;
                    if value.is_empty() {
                        Err(Error::InputStringError("Missing parameter value"))
                    } else {
                        Ok((name, Some(value), rest))
                    }
                } else {
                    // no value parameter
                    Ok((name, None, rest))
                }
            }
        } else {
            // Swallow non parameters parts
            Self::unquote_next(s).map(|(_, rest)| ("", None, rest))
        }
    }

    pub(crate) fn parse(s: &str) -> Result<()> {
        Ok(()) 
    
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projstring_unquote() {
        let (s, r) = Parser::unquote_next("foo").unwrap();
        assert_eq!((s, r), ("foo", ""));

        let s = r#"foo " foobar" foo""bar  "foo ""bar" "baz "#;
        let (s, r) = Parser::unquote_next(s).unwrap();
        assert_eq!((s, r), ("foo", r#"" foobar" foo""bar  "foo ""bar" "baz "#));
        let (s, r) = Parser::unquote_next(r).unwrap();
        assert_eq!((s, r), (" foobar", r#" foo""bar  "foo ""bar" "baz "#));
        let (s, r) = Parser::unquote_next(r).unwrap();
        assert_eq!((s, r), (r#"foo""bar"#, r#" "foo ""bar" "baz "#));
        let (s, r) = Parser::unquote_next(r).unwrap();
        assert_eq!((s, r), (r#"foo ""bar"#, r#" "baz "#));

        assert!(Parser::unquote_next(r).is_err());
    }

    #[test]
    fn projstring_invalid_parameter_name() {
        let s = "+pro@j=geocent";
        assert!(Parser::token(s).is_err());
    }

    #[test]
    fn projstring_token() {
        let s = "+proj=geocent +datum=WGS84 +no_defs";
        let r = Parser::token(s).unwrap();
        assert_eq!(r, ("proj", Some("geocent"), "+datum=WGS84 +no_defs"));
        let r = Parser::token(r.2).unwrap();
        assert_eq!(r, ("datum", Some("WGS84"), "+no_defs"));
        let r = Parser::token(r.2).unwrap();
        assert_eq!(r, ("no_defs", None, ""));
    }
}
