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
use crate::errors::{Error, Result};
use crate::parameters::{ParamList, Parameter};

pub fn parse(s: &str) -> Result<ParamList<'_>> {
    tokenizer::tokens(s)
        .map(|r| match r {
            Ok((name, value, _)) => Ok(Parameter { name, value }),
            Err(err) => Err(err),
        })
        .collect()
}

mod tokenizer {
    use super::*;
    use std::ops::ControlFlow;

    /// Parse parameter name as valid identifier
    ///
    /// Return error if the string is not a valid
    /// identifier: i.e not [0-9a-zA-Z_]+
    pub(super) fn parse_identifier(s: &str) -> Result<(&str, &str)> {
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
    pub(super) fn unquote_next(s: &str) -> Result<(&str, &str)> {
        let s = s.trim_start();
        if s.starts_with('\"') {
            // Check if string part is terminated by a quote,
            // if not, continue with the next part.
            // Inner quotes not separated by a whitespace is left
            // as part of the token.
            let s = unsafe { s.split_once('\"').unwrap_unchecked().1 };
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
                ControlFlow::Continue(_) => {
                    Err(Error::InputStringError("Unterminated quoted string"))
                }
            }
        } else {
            Ok(s.split_once(|c: char| c.is_whitespace()).unwrap_or((s, "")))
        }
    }

    /// Returns the first token from the input str
    pub(super) fn token(s: &str) -> Result<(&str, Option<&str>, &str)> {
        let s = s.trim_start();
        if s.is_empty() {
            Ok(("", None, ""))
        } else if s.starts_with('+') {
            let (_, rest) = unsafe { s.split_once('+').unwrap_unchecked() }; // Swallow '+'

            let (name, rest) = parse_identifier(rest)?;
            if name.is_empty() {
                Err(Error::InputStringError("Empty parameter name"))
            } else {
                let rest = rest.trim_start();
                if rest.starts_with('=') {
                    let (value, rest) = unquote_next(rest.split_once('=').unwrap().1)?;
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
            unquote_next(s).map(|(_, rest)| ("", None, rest))
        }
    }

    /// Generate an iterator from parsing results
    pub(super) fn tokens(s: &str) -> impl Iterator<Item = Result<(&str, Option<&str>, &str)>> {
        std::iter::successors(
            Some(token(s)),
            |prev: &Result<(&str, Option<&str>, &str)>| match prev {
                Err(_) => None,
                Ok((_, _, s)) => {
                    if s.is_empty() {
                        None
                    } else {
                        Some(tokenizer::token(s))
                    }
                }
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test tokenizer

    #[test]
    fn projstring_unquote() {
        let (s, r) = tokenizer::unquote_next("foo").unwrap();
        assert_eq!((s, r), ("foo", ""));

        let s = r#"foo " foobar" foo""bar  "foo ""bar" "baz "#;
        let (s, r) = tokenizer::unquote_next(s).unwrap();
        assert_eq!((s, r), ("foo", r#"" foobar" foo""bar  "foo ""bar" "baz "#));
        let (s, r) = tokenizer::unquote_next(r).unwrap();
        assert_eq!((s, r), (" foobar", r#" foo""bar  "foo ""bar" "baz "#));
        let (s, r) = tokenizer::unquote_next(r).unwrap();
        assert_eq!((s, r), (r#"foo""bar"#, r#" "foo ""bar" "baz "#));
        let (s, r) = tokenizer::unquote_next(r).unwrap();
        assert_eq!((s, r), (r#"foo ""bar"#, r#" "baz "#));

        assert!(tokenizer::unquote_next(r).is_err());
    }

    #[test]
    fn projstring_invalid_parameter_name() {
        let s = "+pro@j=geocent";
        assert!(tokenizer::token(s).is_err());
    }

    #[test]
    fn projstring_token() {
        let s = "+proj=geocent +datum=WGS84 +no_defs";
        let r = tokenizer::token(s).unwrap();
        assert_eq!(r, ("proj", Some("geocent"), "+datum=WGS84 +no_defs"));
        let r = tokenizer::token(r.2).unwrap();
        assert_eq!(r, ("datum", Some("WGS84"), "+no_defs"));
        let r = tokenizer::token(r.2).unwrap();
        assert_eq!(r, ("no_defs", None, ""));
    }

    #[test]
    fn projstring_collect_tokens() {
        // Check valid projstring
        let s = "+proj=geocent +datum=WGS84 +no_defs";
        let tokens: Result<Vec<_>> = tokenizer::tokens(s).collect();
        assert!(tokens.is_ok());
        assert_eq!(tokens.unwrap().len(), 3);

        // Check invalid parameters
        let s = "+pro@j=geocent";
        let tokens: Result<Vec<_>> = tokenizer::tokens(s).collect();
        assert!(tokens.is_err())
    }
}
