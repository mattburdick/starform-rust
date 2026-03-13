/// Small profiling harness that generates systems for the HigherMultiplicity
/// seeds (1, 14, 15) repeatedly so a sampling profiler can capture hotspots.
use starform_rust::{
    accretion_parameters::{set_all_accretion_parameters, AccretionParameters},
    consts::{DUST_DENSITY_COEFF, K},
    generate_star_system, loglevel_mutex,
    random::set_rng_seed,
};

fn main() {
    // Silence logging
    *loglevel_mutex().lock().unwrap() = 0;

    // Use defaults matching the baseline test
    set_all_accretion_parameters(AccretionParameters {
        dust_density_coefficient: DUST_DENSITY_COEFF,
        percent_dust_in_cloud: K,
        ..AccretionParameters::default()
    });

    let seeds: &[u64] = &[1, 14, 15, 3, 4, 8, 2, 5, 6];
    let iterations = 10;

    for iteration in 0..iterations {
        for &seed in seeds {
            set_rng_seed(seed);
            let system = generate_star_system(0, "".to_owned());
            // Prevent optimizing away
            std::hint::black_box(&system);
        }
        eprintln!("iteration {}/{iterations}", iteration + 1);
    }
}
