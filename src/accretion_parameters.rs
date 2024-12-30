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
    /// println!("Percent of cloud as dust: {}", params.percent_dust_in_cloud);
    /// ```
    pub static ref ACCRETION_PARAMETERS: Mutex<AccretionParameters> = Mutex::new(AccretionParameters {
        dust_density_coefficient: DUST_DENSITY_COEFF,
        percent_dust_in_cloud: ((1.0 / K)*100.0).min(100.0),
    });
}

/// A struct to hold accretion parameters, including dust density coefficient and alpha.
/// These parameters are used in density calculations and other related computations.
pub struct AccretionParameters {
    /// "Parameter A" in Dole's paper
    pub dust_density_coefficient: f64,

    /// 1/"Parameter K" in Dole's paper. K is the ratio of gas to dust, typically in the range 50 to 100, so percent_dust_in_cloud
    /// should typically range from 1% to 2%.
    pub percent_dust_in_cloud: f64,
}

/// Retrieves the current accretion parameters by locking the global mutex.
///
/// # Returns
/// An instance of `AccretionParameters` representing the current values.
pub fn get_accretion_parameters() -> AccretionParameters {
    let params = ACCRETION_PARAMETERS.lock().unwrap();
    AccretionParameters {
        dust_density_coefficient: params.dust_density_coefficient,
        percent_dust_in_cloud: params.percent_dust_in_cloud,
    }
}

/// Updates the accretion parameters by locking the global mutex and modifying the fields.
///
/// # Arguments
/// - `new_dust_density_coeff`: The new dust density coefficient.
/// - `percent_dust_in_cloud`: The new percent of the cloud that is dust.
pub fn set_accretion_parameters(new_dust_density_coeff: f64, new_percent_dust_in_cloud: f64) {
    let mut params = ACCRETION_PARAMETERS.lock().unwrap();

    // Assume we only want to update a setting if the new setting is non-zero
    if new_dust_density_coeff != 0.0 {
        params.dust_density_coefficient = new_dust_density_coeff;
    }

    if new_percent_dust_in_cloud != 0.0 {
        params.percent_dust_in_cloud = new_percent_dust_in_cloud;
    }
}
