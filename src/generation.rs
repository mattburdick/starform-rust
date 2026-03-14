use crate::{
    accretion_parameters::AccretionParameters,
    star::{LuminosityClass, SpectralClass},
};

pub const MAX_GENERATION_STARS: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SystemGenerationPreset {
    Custom,
    RockyRich,
    SolarLike,
    GiantRich,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostStarPreset {
    Custom,
    RandomObserved,
    CoolDwarf,
    SolarAnalog,
    WarmBright,
    GiantHost,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiplicityPreset {
    Custom,
    Observed,
    SingleHeavy,
    BinaryHeavy,
    HigherMultiplicity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StarCountMode {
    Sampled,
    Exact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeparationMode {
    Sampled,
    Exact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StarSelectionMode {
    RandomObserved,
    Explicit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConfiguredStar {
    pub selection_mode: StarSelectionMode,
    pub spectral_class: SpectralClass,
    pub spectral_number: i32,
    pub luminosity_class: LuminosityClass,
}

impl ConfiguredStar {
    pub fn explicit(spectral_class: SpectralClass, spectral_number: i32, luminosity_class: LuminosityClass) -> Self {
        Self {
            selection_mode: StarSelectionMode::Explicit,
            spectral_class,
            spectral_number,
            luminosity_class,
        }
    }
}

impl Default for ConfiguredStar {
    fn default() -> Self {
        Self {
            selection_mode: StarSelectionMode::RandomObserved,
            spectral_class: SpectralClass::G,
            spectral_number: 2,
            luminosity_class: LuminosityClass::MainSequence,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SystemGenerationConfig {
    pub system_preset: SystemGenerationPreset,
    pub host_star_preset: HostStarPreset,
    pub multiplicity_preset: MultiplicityPreset,
    pub star_count_mode: StarCountMode,
    pub star_count: usize,
    pub separation_mode: SeparationMode,
    pub binary_separation_au: f64,
    pub stars: [ConfiguredStar; MAX_GENERATION_STARS],
    pub accretion: AccretionParameters,
}

impl Default for SystemGenerationConfig {
    fn default() -> Self {
        Self {
            system_preset: SystemGenerationPreset::Custom,
            host_star_preset: HostStarPreset::RandomObserved,
            multiplicity_preset: MultiplicityPreset::Observed,
            star_count_mode: StarCountMode::Sampled,
            star_count: 1,
            separation_mode: SeparationMode::Sampled,
            binary_separation_au: 24.0,
            stars: [ConfiguredStar::default(); MAX_GENERATION_STARS],
            accretion: AccretionParameters::default(),
        }
    }
}

impl SystemGenerationConfig {
    pub fn rocky_rich() -> Self {
        let mut config = Self {
            system_preset: SystemGenerationPreset::RockyRich,
            host_star_preset: HostStarPreset::CoolDwarf,
            multiplicity_preset: MultiplicityPreset::SingleHeavy,
            star_count_mode: StarCountMode::Exact,
            star_count: 1,
            separation_mode: SeparationMode::Sampled,
            binary_separation_au: 24.0,
            stars: [ConfiguredStar::default(); MAX_GENERATION_STARS],
            accretion: AccretionParameters {
                dust_density_coefficient: 7.8e-8,
                percent_dust_in_cloud: 2.4,
                disk_mass_multiplier: 0.85,
                metallicity: 1.30,
                gas_disk_lifetime_multiplier: 0.55,
                giant_core_efficiency: 0.65,
            },
        };
        config.stars[0] = ConfiguredStar::explicit(SpectralClass::K, 2, LuminosityClass::MainSequence);
        config
    }

    pub fn solar_like() -> Self {
        let mut config = Self {
            system_preset: SystemGenerationPreset::SolarLike,
            host_star_preset: HostStarPreset::SolarAnalog,
            multiplicity_preset: MultiplicityPreset::Observed,
            star_count_mode: StarCountMode::Exact,
            star_count: 1,
            separation_mode: SeparationMode::Sampled,
            binary_separation_au: 24.0,
            stars: [ConfiguredStar::default(); MAX_GENERATION_STARS],
            accretion: AccretionParameters::default(),
        };
        config.stars[0] = ConfiguredStar::explicit(SpectralClass::G, 2, LuminosityClass::MainSequence);
        config
    }

    pub fn giant_rich() -> Self {
        let mut config = Self {
            system_preset: SystemGenerationPreset::GiantRich,
            host_star_preset: HostStarPreset::WarmBright,
            multiplicity_preset: MultiplicityPreset::SingleHeavy,
            star_count_mode: StarCountMode::Exact,
            star_count: 1,
            separation_mode: SeparationMode::Sampled,
            binary_separation_au: 24.0,
            stars: [ConfiguredStar::default(); MAX_GENERATION_STARS],
            accretion: AccretionParameters {
                dust_density_coefficient: 1.1e-7,
                percent_dust_in_cloud: 2.0,
                disk_mass_multiplier: 1.70,
                metallicity: 1.20,
                gas_disk_lifetime_multiplier: 1.85,
                giant_core_efficiency: 1.75,
            },
        };
        config.stars[0] = ConfiguredStar::explicit(SpectralClass::F, 2, LuminosityClass::MainSequence);
        config
    }

    pub fn apply_system_preset(&mut self, preset: SystemGenerationPreset) {
        *self = match preset {
            SystemGenerationPreset::Custom => Self::default(),
            SystemGenerationPreset::RockyRich => Self::rocky_rich(),
            SystemGenerationPreset::SolarLike => Self::solar_like(),
            SystemGenerationPreset::GiantRich => Self::giant_rich(),
        };
    }

    pub fn apply_host_star_preset(&mut self, preset: HostStarPreset) {
        self.system_preset = SystemGenerationPreset::Custom;
        self.host_star_preset = preset;
        self.stars[0] = match preset {
            HostStarPreset::Custom => self.stars[0],
            HostStarPreset::RandomObserved => ConfiguredStar {
                selection_mode: StarSelectionMode::RandomObserved,
                ..self.stars[0]
            },
            HostStarPreset::CoolDwarf => ConfiguredStar::explicit(SpectralClass::M, 3, LuminosityClass::MainSequence),
            HostStarPreset::SolarAnalog => ConfiguredStar::explicit(SpectralClass::G, 2, LuminosityClass::MainSequence),
            HostStarPreset::WarmBright => ConfiguredStar::explicit(SpectralClass::F, 2, LuminosityClass::MainSequence),
            HostStarPreset::GiantHost => ConfiguredStar::explicit(SpectralClass::K, 0, LuminosityClass::Giant),
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn giant_rich_preset_boosts_giant_friendly_knobs() {
        let config = SystemGenerationConfig::giant_rich();
        assert_eq!(config.system_preset, SystemGenerationPreset::GiantRich);
        assert_eq!(config.star_count, 1);
        assert!(config.accretion.disk_mass_multiplier > 1.0);
        assert!(config.accretion.gas_disk_lifetime_multiplier > 1.0);
        assert!(config.accretion.giant_core_efficiency > 1.0);
    }

    #[test]
    fn host_star_preset_sets_primary_override() {
        let mut config = SystemGenerationConfig::default();
        config.apply_host_star_preset(HostStarPreset::SolarAnalog);
        assert_eq!(config.host_star_preset, HostStarPreset::SolarAnalog);
        assert_eq!(config.stars[0].selection_mode, StarSelectionMode::Explicit);
        assert_eq!(config.stars[0].spectral_class, SpectralClass::G);
        assert_eq!(config.stars[0].spectral_number, 2);
    }
}
