//!
//! Display benchmarks for computing forward and inverse projections
//!
//! Compute benchmarks the same way as the `bench_proj_trans` PROJ
//! test utility
//!
use clap::{ArgAction, Parser};
use proj4rs::{
    errors::{Error, Result},
    proj, transform,
};

use rand::prelude::*;

use std::io::{self, BufRead};
use std::time::Instant;

#[derive(Parser)]
#[command(author, version="0.1", about="Bench projections", long_about = None)]
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
    #[arg(short, long, default_value_t = 5_000_000)]
    loops: u32,
    #[arg(long, default_value_t = 0.0)]
    noise_x: f64,
    #[arg(long, default_value_t = 0.0)]
    noise_y: f64,
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

    let loops = args.loops;
    let noise_x = args.noise_x;
    let noise_y = args.noise_y;

    let src = proj::Proj::from_user_string(srcdef)?;
    let dst = proj::Proj::from_user_string(dstdef)?;

    let stdin = io::stdin().lock();

    fn from_parse_err(err: std::num::ParseFloatError) -> Error {
        Error::ParameterValueError(format!("{err:?}"))
    }

    for line in stdin.lines() {
        let line = line.unwrap();
        let inputs = line.as_str().split_whitespace().collect::<Vec<_>>();
        if inputs.len() < 2 || inputs.len() > 3 {
            eprintln!("Expecting: '<x> ,<y> [,<z>]' found: {}", line.as_str());
            std::process::exit(1);
        }

        let mut x: f64 = inputs[0].parse().map_err(from_parse_err)?;
        let mut y: f64 = inputs[1].parse().map_err(from_parse_err)?;
        let z: f64 = if inputs.len() > 2 {
            inputs[2].parse().map_err(from_parse_err)?
        } else {
            0.
        };

        if src.is_latlong() {
            x = x.to_radians();
            y = y.to_radians();
        }

        let mut point = (x, y, z);

        transform::transform(&src, &dst, &mut point)?;
        println!("{} {} -> {:.16} {:.16}", x, y, point.0, point.1);

        let mut rng = rand::rng();

        // Time noise generation
        let start = Instant::now();
        for _ in 0..loops {
            if noise_x > 0. {
                point.0 = x + noise_x * (2. * rng.random::<f64>() - 1.);
            }
            if noise_y > 0. {
                point.1 = y + noise_y * (2. * rng.random::<f64>() - 1.);
            }
        }
        let noise_elapsed = start.elapsed();

        let start = Instant::now();
        for _ in 0..loops {
            if noise_x > 0. {
                point.0 = x + noise_x * (2. * rng.random::<f64>() - 1.);
            }
            if noise_y > 0. {
                point.1 = y + noise_y * (2. * rng.random::<f64>() - 1.);
            }
            transform::transform(&src, &dst, &mut point)?;
        }

        let elapsed = start.elapsed();

        println!("Duration: {} ms", (elapsed - noise_elapsed).as_millis());
        println!(
            "Throughput: {:.2} million coordinates/s",
            1e-3 * (loops as f64) / (elapsed - noise_elapsed).as_millis() as f64
        );
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
