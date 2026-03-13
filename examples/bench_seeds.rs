/// Benchmark harness that times system generation for each seed bucket.
/// Reports per-bucket mean times in milliseconds.
use starform_rust::{
    accretion_parameters::{set_all_accretion_parameters, AccretionParameters},
    consts::{DUST_DENSITY_COEFF, K},
    generate_star_system, loglevel_mutex,
    random::set_rng_seed,
};
use std::time::Instant;

fn measure_seed(seed: u64, runs: usize) -> f64 {
    // Warmup
    set_rng_seed(seed);
    let _ = std::hint::black_box(generate_star_system(0, "".to_owned()));

    let mut times = Vec::with_capacity(runs);
    for _ in 0..runs {
        set_rng_seed(seed);
        let start = Instant::now();
        let system = generate_star_system(0, "".to_owned());
        let elapsed = start.elapsed();
        std::hint::black_box(&system);
        times.push(elapsed.as_secs_f64() * 1000.0);
    }
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // Trimmed mean: drop min and max if we have enough runs
    if times.len() > 2 {
        let trimmed = &times[1..times.len() - 1];
        trimmed.iter().sum::<f64>() / trimmed.len() as f64
    } else {
        times.iter().sum::<f64>() / times.len() as f64
    }
}

fn main() {
    *loglevel_mutex().lock().unwrap() = 0;

    set_all_accretion_parameters(AccretionParameters {
        dust_density_coefficient: DUST_DENSITY_COEFF,
        percent_dust_in_cloud: K,
        ..AccretionParameters::default()
    });

    let buckets: &[(&str, &[u64])] = &[
        ("SingleStar", &[2, 5, 6]),
        ("Binary", &[3, 4, 8]),
        ("HigherMultiplicity", &[1, 14, 15]),
    ];

    let runs = 5;

    eprintln!("Running {} timed iterations per seed (+ warmup)...\n", runs);

    for (name, seeds) in buckets {
        let mut bucket_times = Vec::new();
        for &seed in *seeds {
            let mean_ms = measure_seed(seed, runs);
            eprintln!("  seed {seed:>2}: {mean_ms:>8.2} ms");
            bucket_times.push(mean_ms);
        }
        let bucket_mean = bucket_times.iter().sum::<f64>() / bucket_times.len() as f64;
        // Budget calculation: round_up(mean * 2.0), floor at 300ms
        let budget = (bucket_mean * 2.0).ceil().max(300.0);
        eprintln!("{name:>20}: mean = {bucket_mean:>8.2} ms, budget = {budget:>8.1} ms");
        eprintln!();
    }
}
