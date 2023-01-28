//!
//! Compare loop vs  try_fold
//!
//!
use proj4rs::adaptors::transform_point_array;
use proj4rs::errors::{Error, Result};
use proj4rs::proj::Proj;

use std::ops::ControlFlow::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn etmerc_transform(itermax: usize) {
    let d = 1.0 / (itermax as f64);

    let mut data: Vec<(f64, f64, f64)> = (1..=itermax)
        .map(|i| {
            (
                (-2.0f64 + (i as f64) * 4.0 * d).to_radians(),
                (-1.0f64 + (i as f64) * 2.0 * d).to_radians(),
                //2.0f64.to_radians(),
                //1.0f64.to_radians(),
                0.,
            )
        })
        .collect();

    let from = Proj::from_proj_string("+proj=latlong +ellps=GRS80").unwrap();
    let to = Proj::from_proj_string("+proj=etmerc +ellps=GRS80").unwrap();

    transform_point_array(&from, &to, data.as_mut_slice()).unwrap();
}

fn criterion_benchmark_proj(c: &mut Criterion) {
    c.bench_function("tmerc forward", |b| {
        b.iter(|| etmerc_transform(black_box(10_000usize)))
    });
}

criterion_group!(benches, criterion_benchmark_proj);
criterion_main!(benches);
