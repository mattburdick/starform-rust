use rand::Rng;

/// Generates a random number within the inclusive range defined by two input numbers.
///
/// This function determines the minimum and maximum values from the provided inputs and then generates
/// a random number within the inclusive range `[min, max]`. It uses the system's default random number generator
/// with a uniform distribution.
///
/// # Parameters
/// - `num1`: One endpoint of the range. This value can be either the minimum or maximum.
/// - `num2`: The other endpoint of the range. This value can be either the minimum or maximum.
///
/// # Returns
/// Returns a random `f64` value that lies between `num1` and `num2`, inclusive of both endpoints.
///
/// # Examples
/// ```
/// let num1 = 10.0;
/// let num2 = 20.0;
/// let random_val = random_number(num1, num2);
/// println!("Random value between {} and {}: {}", num1, num2, random_val);
/// assert!(random_val >= 10.0 && random_val <= 20.0);
/// ```
///
/// # Note
/// The function uses `rand::thread_rng()`, a thread-local random number generator, to produce randomness.
/// This generator is seeded by the system and is safe for casual use but not for cryptographic purposes.
pub fn random_number(num1: f64, num2: f64) -> f64 {
    let max = f64::max(num1, num2);
    let min = f64::min(num1, num2);

    rand::thread_rng().gen_range(min..=max)
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
/// - `percent_variation`: The maximum fractional variation allowed from the base value, expressed as a decimal
///   (e.g., 0.1 for ±10%).
///
/// # Returns
/// Returns a new `f64` that is the original `value` adjusted by a random factor within the specified `percent_variation` range.
///
pub fn about(value: f64, percent_variation: f64) -> f64 {
    let mut rng = rand::thread_rng();
    let random_factor = rng.gen_range(-percent_variation..=percent_variation);
    value + (value * random_factor)
}
