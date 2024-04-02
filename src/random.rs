use rand::Rng;
// use std::cmp;

pub fn random_number(num1: f64, num2: f64) -> f64 {
    let max = f64::max(num1, num2);
    let min = f64::min(num1, num2);

    rand::thread_rng().gen_range(min..=max)
}

/// This function returns a value within a certain variation of the
/// exact value given it in 'value'.
pub fn about(value: f64, variation: f64) -> f64 {
    let mut rng = rand::thread_rng();
    let random_factor = rng.gen_range(-variation..=variation);
    value + (value * random_factor)
}
