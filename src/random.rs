use rand::{rngs::StdRng, Rng, SeedableRng};
use std::{ops::RangeInclusive, sync::Mutex};

// A globally accessible, thread-safe random number generator (RNG).
//
// This RNG is initialized once using a fixed seed, ensuring deterministic random number generation across runs. The RNG
// is wrapped in a `Mutex` to allow safe, synchronized access across threads.
lazy_static! {
    /// A globally shared, thread-safe `StdRng` initialized with a fixed seed.
    static ref DETERMINISTIC_RNG: Mutex<StdRng> = Mutex::new(StdRng::seed_from_u64(rand::thread_rng().gen()));
}

/// Generate a new random seed and reset the global RNG with it.
///
/// # Returns
/// - The new random seed (`u64`) that was generated.
///
/// # Example
/// ```
/// let new_seed = reset_rng();
/// println!("RNG reset with new seed: {}", new_seed);
/// ```
pub fn set_rng_seed(seed: u64) -> u64 {
    // Generate a new random seed using the default thread_rng()
    let new_seed: u64 = if seed == 0 { rand::thread_rng().gen() } else { seed };

    // let log_level = *crate::get_log_level!();
    // crate::log!(log_level, 1, "Seeding with {}", new_seed);

    // Lock the RNG and replace it with a new StdRng seeded with the new seed
    let mut rng = DETERMINISTIC_RNG.lock().unwrap();
    *rng = StdRng::seed_from_u64(new_seed);

    new_seed
}

/// Generates a random number of a specified numeric type in the given range `[min, max]`
/// using the global deterministic RNG.
///
/// # Arguments
/// - `range`: A range that implements `std::ops::RangeInclusive<T>`.
///
/// # Type Parameters
/// - `T`: A numeric type that implements `rand::distributions::uniform::SampleUniform`.
///
/// # Returns
/// - A random number of type `T` within the specified range.
pub fn get_random_number<T>(range: RangeInclusive<T>) -> T
where
    T: rand::distributions::uniform::SampleUniform + PartialOrd,
{
    DETERMINISTIC_RNG.lock().unwrap().gen_range(range)
}

/// Calculates a value within a specific range around a given value, introducing random variation.
///
/// This function generates a new value by adding a random percentage of the original `value` to itself. The percentage
/// is randomly chosen within the range of `-percent_variation` to `percent_variation`, where `percent_variation` is expressed as a decimal
/// percentage of the original value. For example, a `percent_variation` of 0.1 allows the function to return a value between
/// 90% and 110% of the original value.
///
/// # Parameters
/// - `value`: The base value from which the variation is calculated.
/// - `percent_variation`: The maximum fractional v÷≥ariation allowed from the base value, expressed as a decimal
///   (e.g., 0.1 for ±10%).
///
/// # Returns
/// Returns a new `f64` that is the original `value` adjusted by a random factor within the specified `percent_variation` range.
///
pub fn about(value: f64, percent_variation: f64) -> f64 {
    let mut rng = DETERMINISTIC_RNG.lock().unwrap();
    let random_factor = rng.gen_range(-percent_variation..=percent_variation);
    value + (value * random_factor)
}
