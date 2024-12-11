use std::sync::Mutex;

use crate::consts::{DUST_DENSITY_COEFF, K};

lazy_static! {
    /// A global, thread-safe set of accretion parameters.
    ///
    /// The `ACCRETION_PARAMETERS` struct stores critical values related to density calculations and
    /// is guarded by a `Mutex` to ensure safe concurrent access.
    ///
    /// Access the fields by locking the mutex:
    /// ```rust
    /// let params = ACCRETION_PARAMETERS.lock().unwrap();
    /// println!("Dust Density Coeff: {}", params.dust_density_coeff);
    /// println!("Ratio of gas to dust: {}", params.ratio_of_gas_to_dust);
    /// ```
    pub static ref ACCRETION_PARAMETERS: Mutex<AccretionParameters> = Mutex::new(AccretionParameters {
        dust_density_coefficient: DUST_DENSITY_COEFF,
        ratio_of_gas_to_dust: K,
    });
}

/// A struct to hold accretion parameters, including dust density coefficient and alpha.
/// These parameters are used in density calculations and other related computations.
pub struct AccretionParameters {
    /// "Parameter A" in Dole's paper
    pub dust_density_coefficient: f64,

    /// "Parameter K" in Dole's paper. The ratio of gas to dust, typically in the range 50 to 100 (e.g. 1% to 2% of the cloud mass is dust)
    pub ratio_of_gas_to_dust: f64,
}

/// Retrieves the current accretion parameters by locking the global mutex.
///
/// # Returns
/// An instance of `AccretionParameters` representing the current values.
pub fn get_accretion_parameters() -> AccretionParameters {
    let params = ACCRETION_PARAMETERS.lock().unwrap();
    AccretionParameters {
        dust_density_coefficient: params.dust_density_coefficient,
        ratio_of_gas_to_dust: params.ratio_of_gas_to_dust,
    }
}

/// Updates the accretion parameters by locking the global mutex and modifying the fields.
///
/// # Arguments
/// - `new_dust_density_coeff`: The new dust density coefficient.
/// - `new_ratio_of_gas_to_dust`: The new ratio of gas to dust.
pub fn set_accretion_parameters(new_dust_density_coeff: f64, new_ratio_of_gas_to_dust: f64) {
    let mut params = ACCRETION_PARAMETERS.lock().unwrap();

    // Assume we only want to update a setting if the new setting is non-zero
    if new_dust_density_coeff != 0.0 {
        params.dust_density_coefficient = new_dust_density_coeff;
    }

    if new_ratio_of_gas_to_dust != 0.0 {
        params.ratio_of_gas_to_dust = new_ratio_of_gas_to_dust;
    }
}
