//!
//! Display debugging projection infos
//!
use proj4rs::{
    errors::{Error, Result},
    proj,
};
use std::env;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 {
        println!("Usage: projdbg <projstr>");
        return Err(Error::InvalidParameterValue("Missing proj string"));
    }

    let projstr = args[1..].join(" ");
    let projection = proj::Proj::from_user_string(&projstr)?;

    println!("{projection:#?}");
    Ok(())
}
