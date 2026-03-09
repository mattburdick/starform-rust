#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CondensationThresholds {
    pub refractory_metal_temp_k: f64,
    pub silicate_rock_temp_k: f64,
    pub water_ice_temp_k: f64,
    pub volatile_ice_temp_k: f64,
}

impl CondensationThresholds {
    pub fn is_descending(&self) -> bool {
        self.refractory_metal_temp_k > self.silicate_rock_temp_k
            && self.silicate_rock_temp_k > self.water_ice_temp_k
            && self.water_ice_temp_k > self.volatile_ice_temp_k
    }
}

impl Default for CondensationThresholds {
    fn default() -> Self {
        Self {
            refractory_metal_temp_k: 1_350.0,
            silicate_rock_temp_k: 950.0,
            water_ice_temp_k: 180.0,
            volatile_ice_temp_k: 70.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegionThermalProfile {
    pub host_luminosity_solar: f64,
    pub reference_temperature_k: f64,
    pub radial_exponent: f64,
    pub transition_width_k: f64,
}

impl RegionThermalProfile {
    pub fn from_host_luminosity(host_luminosity_solar: f64) -> Self {
        let luminosity = host_luminosity_solar.max(0.01);
        Self {
            host_luminosity_solar,
            reference_temperature_k: 278.0 * luminosity.powf(0.25),
            radial_exponent: 0.5,
            transition_width_k: 25.0,
        }
    }

    pub fn temperature_at_au(&self, orbital_radius_au: f64) -> f64 {
        let radius = orbital_radius_au.max(0.01);
        self.reference_temperature_k / radius.powf(self.radial_exponent)
    }

    pub fn sample_planet_formation_inputs(
        &self,
        orbital_radius_au: f64,
        thresholds: CondensationThresholds,
    ) -> PlanetFormationInputs {
        debug_assert!(thresholds.is_descending());
        let temperature_k = self.temperature_at_au(orbital_radius_au);
        let condensation_fractions =
            CondensationFractions::from_temperature(temperature_k, thresholds, self.transition_width_k);

        PlanetFormationInputs {
            temperature_k,
            condensation_fractions,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CondensationFractions {
    pub refractory_metal: f64,
    pub silicate_rock: f64,
    pub water_ice: f64,
    pub volatile_ices: f64,
    pub gas: f64,
}

impl CondensationFractions {
    pub fn from_temperature(temperature_k: f64, thresholds: CondensationThresholds, transition_width_k: f64) -> Self {
        let refractory_metal = weighted_fraction(
            temperature_k,
            thresholds.refractory_metal_temp_k,
            transition_width_k,
            0.15,
        );
        let silicate_rock = weighted_fraction(temperature_k, thresholds.silicate_rock_temp_k, transition_width_k, 0.35);
        let water_ice = weighted_fraction(temperature_k, thresholds.water_ice_temp_k, transition_width_k, 0.30);
        let volatile_ices = weighted_fraction(temperature_k, thresholds.volatile_ice_temp_k, transition_width_k, 0.20);
        let gas = (1.0 - (refractory_metal + silicate_rock + water_ice + volatile_ices)).clamp(0.0, 1.0);

        Self {
            refractory_metal,
            silicate_rock,
            water_ice,
            volatile_ices,
            gas,
        }
    }

    pub fn total(&self) -> f64 {
        self.refractory_metal + self.silicate_rock + self.water_ice + self.volatile_ices + self.gas
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanetFormationInputs {
    pub temperature_k: f64,
    pub condensation_fractions: CondensationFractions,
}

fn weighted_fraction(temperature_k: f64, threshold_k: f64, transition_width_k: f64, max_fraction: f64) -> f64 {
    max_fraction * smooth_activation_below_threshold(temperature_k, threshold_k, transition_width_k)
}

fn smooth_activation_below_threshold(temperature_k: f64, threshold_k: f64, transition_width_k: f64) -> f64 {
    if transition_width_k <= 0.0 {
        return if temperature_k <= threshold_k { 1.0 } else { 0.0 };
    }

    let half_width = transition_width_k / 2.0;
    let delta = threshold_k - temperature_k;

    if delta <= -half_width {
        0.0
    } else if delta >= half_width {
        1.0
    } else {
        (0.5 + (delta / transition_width_k)).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1.0e-9;

    #[test]
    fn default_thresholds_are_descending() {
        assert!(CondensationThresholds::default().is_descending());
    }

    #[test]
    fn temperature_drops_with_orbital_distance() {
        let profile = RegionThermalProfile::from_host_luminosity(1.0);
        assert!(profile.temperature_at_au(0.5) > profile.temperature_at_au(2.0));
    }

    #[test]
    fn hotter_inner_disk_keeps_more_material_in_gas() {
        let profile = RegionThermalProfile::from_host_luminosity(1.0);
        let thresholds = CondensationThresholds::default();
        let inner = profile.sample_planet_formation_inputs(0.25, thresholds);
        let outer = profile.sample_planet_formation_inputs(8.0, thresholds);

        assert!(inner.condensation_fractions.gas > outer.condensation_fractions.gas);
    }

    #[test]
    fn cooler_outer_disk_accumulates_more_ices() {
        let profile = RegionThermalProfile::from_host_luminosity(1.0);
        let thresholds = CondensationThresholds::default();
        let inner = profile.sample_planet_formation_inputs(0.5, thresholds);
        let outer = profile.sample_planet_formation_inputs(20.0, thresholds);

        assert!(outer.condensation_fractions.water_ice > inner.condensation_fractions.water_ice);
        assert!(outer.condensation_fractions.volatile_ices >= inner.condensation_fractions.volatile_ices);
        assert!((outer.condensation_fractions.total() - 1.0).abs() < EPSILON);
    }

    #[test]
    fn inner_and_outer_sampling_lock_in_distinct_composition_regimes() {
        let profile = RegionThermalProfile::from_host_luminosity(1.0);
        let thresholds = CondensationThresholds::default();
        let inner = profile.sample_planet_formation_inputs(0.01, thresholds);
        let outer = profile.sample_planet_formation_inputs(30.0, thresholds);

        assert!(inner.temperature_k > thresholds.refractory_metal_temp_k);
        assert!(outer.temperature_k < thresholds.volatile_ice_temp_k);

        assert!((inner.condensation_fractions.gas - 1.0).abs() < EPSILON);
        assert!(inner.condensation_fractions.refractory_metal.abs() < EPSILON);
        assert!(inner.condensation_fractions.silicate_rock.abs() < EPSILON);
        assert!(inner.condensation_fractions.water_ice.abs() < EPSILON);
        assert!(inner.condensation_fractions.volatile_ices.abs() < EPSILON);

        assert!(outer.condensation_fractions.gas.abs() < EPSILON);
        assert!((outer.condensation_fractions.refractory_metal - 0.15).abs() < EPSILON);
        assert!((outer.condensation_fractions.silicate_rock - 0.35).abs() < EPSILON);
        assert!((outer.condensation_fractions.water_ice - 0.30).abs() < EPSILON);
        assert!((outer.condensation_fractions.volatile_ices - 0.20).abs() < EPSILON);
    }

    #[test]
    fn repeated_sampling_is_deterministic() {
        let profile = RegionThermalProfile::from_host_luminosity(0.8);
        let thresholds = CondensationThresholds::default();
        let first = profile.sample_planet_formation_inputs(3.2, thresholds);
        let second = profile.sample_planet_formation_inputs(3.2, thresholds);

        assert_eq!(first, second);
    }

    #[test]
    fn transition_width_softens_threshold_edges() {
        let thresholds = CondensationThresholds::default();
        let fractions = CondensationFractions::from_temperature(thresholds.water_ice_temp_k, thresholds, 40.0);

        assert!(fractions.water_ice > 0.0);
        assert!(fractions.water_ice < 0.30);
    }
}
