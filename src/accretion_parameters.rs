use std::sync::Mutex;

use crate::consts::{DUST_DENSITY_COEFF, K};

fn default_percent_dust_in_cloud() -> f64 {
    ((1.0 / K) * 100.0).min(100.0)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AccretionParameters {
    /// "Parameter A" in Dole's paper
    pub dust_density_coefficient: f64,

    /// 1/"Parameter K" in Dole's paper. K is the ratio of gas to dust, typically in the range 50 to 100, so percent_dust_in_cloud
    /// should typically range from 1% to 2%.
    pub percent_dust_in_cloud: f64,

    /// Scales the overall surface density of the disk without changing metallicity.
    pub disk_mass_multiplier: f64,

    /// Scales the solid component of the disk relative to Solar metallicity.
    pub metallicity: f64,

    /// Scales how long nebular gas remains available for envelope accretion.
    pub gas_disk_lifetime_multiplier: f64,

    /// Scales how efficiently large cores can transition into gas-rich planets.
    pub giant_core_efficiency: f64,
}

impl Default for AccretionParameters {
    fn default() -> Self {
        Self {
            dust_density_coefficient: DUST_DENSITY_COEFF,
            percent_dust_in_cloud: default_percent_dust_in_cloud(),
            disk_mass_multiplier: 1.0,
            metallicity: 1.0,
            gas_disk_lifetime_multiplier: 1.0,
            giant_core_efficiency: 1.0,
        }
    }
}

lazy_static! {
    /// A global, thread-safe set of accretion parameters.
    ///
    /// The `ACCRETION_PARAMETERS` struct stores critical values related to density calculations and
    /// is guarded by a `Mutex` to ensure safe concurrent access.
    ///
    /// Access the fields by locking the mutex:
    /// ```rust
    /// use starform_rust::accretion_parameters::ACCRETION_PARAMETERS;
    ///
    /// let params = ACCRETION_PARAMETERS.lock().unwrap();
    /// println!("Dust Density Coeff: {}", params.dust_density_coefficient);
    /// println!("Percent of cloud as dust: {}", params.percent_dust_in_cloud);
    /// ```
    pub static ref ACCRETION_PARAMETERS: Mutex<AccretionParameters> = Mutex::new(AccretionParameters::default());
}

/// Retrieves the current accretion parameters by locking the global mutex.
///
/// # Returns
/// An instance of `AccretionParameters` representing the current values.
pub fn get_accretion_parameters() -> AccretionParameters {
    let params = ACCRETION_PARAMETERS.lock().unwrap();
    *params
}

pub fn set_all_accretion_parameters(new_parameters: AccretionParameters) {
    let mut params = ACCRETION_PARAMETERS.lock().unwrap();
    *params = new_parameters;
}

/// Updates the accretion parameters by locking the global mutex and modifying the fields.
///
/// # Arguments
/// - `new_dust_density_coeff`: The new dust density coefficient.
/// - `percent_dust_in_cloud`: The new percent of the cloud that is dust.
pub fn set_accretion_parameters(new_dust_density_coeff: f64, new_percent_dust_in_cloud: f64) {
    set_all_accretion_parameters(AccretionParameters {
        dust_density_coefficient: new_dust_density_coeff,
        percent_dust_in_cloud: new_percent_dust_in_cloud,
        ..AccretionParameters::default()
    });
}
