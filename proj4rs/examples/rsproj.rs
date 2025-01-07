//!
//! Compute forward and inverse projections
//!
use clap::{ArgAction, Parser};
use proj4rs::{
    errors::{Error, Result},
    proj, transform,
};

use std::io::{self, BufRead};

#[derive(Parser)]
#[command(author, version="0.1", about="Compute projections", long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    /// Destination projection
    #[arg(long, required = true)]
    to: String,
    /// Source projection
    #[arg(long, default_value = "+proj=latlong")]
    from: String,
    /// Perform inverse projection
    #[arg(short, long)]
    inverse: bool,
    /// Increase verbosity
    #[arg(short, long, action = ArgAction::Count)]
    verbose: u8,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    init_logger(args.verbose);

    log::debug!(
        "\nfrom: {}\nto: {}\ninverse: {}",
        args.from,
        args.to,
        args.inverse
    );

    let (srcdef, dstdef): (&str, &str) = if args.inverse {
        (&args.to, &args.from)
    } else {
        (&args.from, &args.to)
    };

    let src = proj::Proj::from_user_string(srcdef)?;
    let dst = proj::Proj::from_user_string(dstdef)?;

    let stdin = io::stdin().lock();

    fn from_parse_err(err: std::num::ParseFloatError) -> Error {
        eprintln!("{err:?}");
        Error::ParameterValueError
    }

    for line in stdin.lines() {
        let line = line.unwrap();
        let inputs = line.as_str().split_whitespace().collect::<Vec<_>>();
        if inputs.len() < 2 || inputs.len() > 3 {
            eprintln!("Expecting: '<x> ,<y> [,<z>]' found: {}", line.as_str());
            std::process::exit(1);
        }

        let x: f64 = inputs[0].parse().map_err(from_parse_err)?;
        let y: f64 = inputs[1].parse().map_err(from_parse_err)?;
        let z: f64 = if inputs.len() > 2 {
            inputs[2].parse().map_err(from_parse_err)?
        } else {
            0.
        };

        let mut point = (x, y, z);

        if src.is_latlong() {
            point.0 = point.0.to_radians();
            point.1 = point.1.to_radians();
        }
        transform::transform(&src, &dst, &mut point)?;
        if dst.is_latlong() {
            point.0 = point.0.to_degrees();
            point.1 = point.1.to_degrees();
        }

        println!("{}  {}", point.0, point.1);
    }
    Ok(())
}

//
// Logger
//
fn init_logger(verbose: u8) {
    use env_logger::Env;
    use log::LevelFilter;

    let mut builder = env_logger::Builder::from_env(Env::default().default_filter_or("info"));

    match verbose {
        1 => builder.filter_level(LevelFilter::Debug),
        _ if verbose > 1 => builder.filter_level(LevelFilter::Trace),
        _ => &mut builder,
    }
    .init();
}
