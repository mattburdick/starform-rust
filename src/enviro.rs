//! Planetary environment routines for determining things like planet radius, density,
//! mass, surface temp, etc.
//!
//! Originally ported from enviro.c (1991), now partially modernized with
//! class-based atmosphere and climate approximations.
//!
//! The atmosphere path in this module is intentionally hybrid: the code still uses the
//! classic Starform/(Fogg 1985)-style helpers for baseline pressure and temperature work, but layers
//! a newer formation-aware atmosphere synthesis step on top.
//!
//! The rationale for that newer step is summarized in `starform-rust/README.md` under
//! `Atmosphere modeling notes`. The papers cited there are the main reference point for the
//! heuristics in this file.

use crate::body::{Body, BodyClimateClass, OrbitalZone};
use crate::consts;
use crate::consts::unused_constants as env;
use crate::random::about;
use crate::star::Star;
use crate::types::MassType;

fn earth_equivalent_insolation_radius_au(star: &Star) -> f64 {
    Star::earth_equivalent_insolation_au(star.luminosity_in_sols)
}

/// Environmental properties for a planet/body.
///
/// These values mirror the fields used in the original C implementation and are
/// intentionally kept together to separate environmental calculations from
/// orbital mechanics.
#[derive(Debug, Clone, Default)]
pub struct EnviroProperties {
    pub a: f64,
    pub radius_in_km: f64,
    pub molec_weight: f64,
    pub mean_molecular_weight: f64,
    pub surf_pressure: f64,
    pub volatile_gas_inventory: f64,
    pub boil_point: f64,
    /// Total water coverage fraction (0–1), liquid and frozen combined.
    /// Use `ice_cover` to determine how much of this is frozen, and
    /// `(hydrosphere - ice_cover).max(0.0)` for the liquid portion.
    pub hydrosphere: f64,
    pub cloud_cover: f64,
    pub ice_cover: f64,
    pub albedo: f64,
    /// Radiative-equilibrium temperature before greenhouse/climate warming.
    pub effective_temp_k: f64,
    /// Additional warming added on top of `effective_temp_k` by the surface model.
    pub greenhouse_rise_k: f64,
    pub surf_temp: f64,
    /// Compact atmosphere label used by tests and future UI/game-facing presentation.
    pub atmosphere_class: AtmosphereClass,
    /// Gas mixture after retention/outgassing heuristics are applied.
    pub atmosphere_components: Vec<AtmosphereGasComponent>,
    /// Dominant gas in `atmosphere_components`, if any.
    pub dominant_gas: Option<AtmosphereGas>,
    /// Derived qualitative traits inferred from the gas mixture.
    pub atmosphere_traits: AtmosphereTraits,
}

/// Coarse atmospheric species tracked by the synthesized gas-mixture model.
///
/// This is intentionally a compact gameplay-oriented set rather than a full photochemical
/// network. The goal is to cover the major regimes called out in the research notes:
/// primordial `H2/He` envelopes, secondary `N2/CO2/H2O` atmospheres, and hot sulfur-bearing
/// corrosive cases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtmosphereGas {
    Hydrogen,
    Helium,
    Methane,
    Ammonia,
    WaterVapor,
    Nitrogen,
    Oxygen,
    CarbonDioxide,
    Argon,
    SulfurDioxide,
}

impl AtmosphereGas {
    pub fn label(self) -> &'static str {
        match self {
            AtmosphereGas::Hydrogen => "H2",
            AtmosphereGas::Helium => "He",
            AtmosphereGas::Methane => "CH4",
            AtmosphereGas::Ammonia => "NH3",
            AtmosphereGas::WaterVapor => "H2O",
            AtmosphereGas::Nitrogen => "N2",
            AtmosphereGas::Oxygen => "O2",
            AtmosphereGas::CarbonDioxide => "CO2",
            AtmosphereGas::Argon => "Ar",
            AtmosphereGas::SulfurDioxide => "SO2",
        }
    }

    pub fn molecular_weight(self) -> f64 {
        match self {
            AtmosphereGas::Hydrogen => env::MOL_HYDROGEN,
            AtmosphereGas::Helium => env::HELIUM,
            AtmosphereGas::Methane => env::METHANE,
            AtmosphereGas::Ammonia => env::AMMONIA,
            AtmosphereGas::WaterVapor => env::WATER_VAPOR,
            AtmosphereGas::Nitrogen => env::MOL_NITROGEN,
            AtmosphereGas::Oxygen => env::MOL_OXYGEN,
            AtmosphereGas::CarbonDioxide => env::CARBON_DIOXIDE,
            AtmosphereGas::Argon => env::ARGON,
            AtmosphereGas::SulfurDioxide => env::SULPH_DIOXIDE,
        }
    }
}

/// One retained atmospheric component after synthesis and species-specific escape filtering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphereGasComponent {
    pub gas: AtmosphereGas,
    pub fraction: f64,
    pub partial_pressure_mb: f64,
    pub retained_inventory: f64,
}

/// Derived qualitative tags inferred from the retained gas mixture.
///
/// These are intentionally broad presentation/gameplay flags, not a full medical or
/// photochemical analysis.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AtmosphereTraits {
    pub breathable: bool,
    pub tainted: bool,
    pub toxic: bool,
    pub oxidizing: bool,
    pub reducing: bool,
    pub corrosive: bool,
    pub steam: bool,
    pub hydrogen_rich: bool,
}

/// Compact atmosphere labels preserved for presentation, testing, and compatibility with the
/// older Starform-style output.
///
/// Internally the richer model works from gas components and derived traits first; this enum is
/// the final summary layer, not the only atmospheric representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AtmosphereClass {
    #[default]
    Airless,
    Trace,
    VeryThin,
    VeryThinTainted,
    Thin,
    ThinTainted,
    Standard,
    StandardTainted,
    Dense,
    DenseTainted,
    Exotic,
    Corrosive,
    Insidious,
    DenseHigh,
    ThinLow,
    Steam,
    H2Rich,
    GasGiant,
    Unusual,
}

/// Legacy-style volatile inventory coefficients.
///
/// The old model selected these entirely from orbital zone. The modernized pipeline still uses
/// the same kind of knobs, but chooses them from climate/formation hints when available.
#[derive(Debug, Clone, Copy)]
struct VolatileRetentionProfile {
    proportion_const: f64,
    retention_divisor: f64,
    greenhouse_candidate: bool,
}

/// Normalized atmosphere facts consumed by `classify_atmosphere`.
///
/// Keeping classification inputs separate makes it clear that `AtmosphereClass` is a compact
/// presentation layer built on top of the richer gas-mixture model.
#[derive(Debug, Clone, Copy)]
struct AtmosphereSignals {
    eq_temp_k: f64,
    smallest_mw_retained: f64,
    mean_molecular_weight: f64,
    volatile_gas_inventory: f64,
    surface_grav_gees: f64,
    pressure_proxy_mb: f64,
    mass_earth: f64,
    gas_fraction: f64,
    water_inventory: f64,
    volatile_inventory_fraction: f64,
    climate_class: Option<BodyClimateClass>,
    greenhouse_effect: bool,
    mass_type: MassType,
    dominant_gas: Option<AtmosphereGas>,
    hydrogen_fraction: f64,
    water_vapor_fraction: f64,
    oxygen_fraction: f64,
    carbon_dioxide_fraction: f64,
    sulfur_dioxide_fraction: f64,
    traits: AtmosphereTraits,
}

/// Detailed atmosphere result before it is mapped back to `AtmosphereClass` and the surface model.
#[derive(Debug, Clone, Default)]
struct AtmosphereModel {
    components: Vec<AtmosphereGasComponent>,
    dominant_gas: Option<AtmosphereGas>,
    traits: AtmosphereTraits,
    volatile_gas_inventory: f64,
    surf_pressure_mb: f64,
    smallest_mw_retained: f64,
    mean_molecular_weight: f64,
    water_inventory: f64,
    volatile_inventory_fraction: f64,
    gas_fraction: f64,
    eq_temp_k: f64,
}

/// Preserves the old zone-based volatile retention coefficients as a fallback.
fn legacy_volatile_retention_profile(zone: OrbitalZone) -> VolatileRetentionProfile {
    match zone {
        OrbitalZone::Zone1 => VolatileRetentionProfile {
            proportion_const: 100000.0,
            retention_divisor: 10.0,
            greenhouse_candidate: true,
        },
        OrbitalZone::Zone2 => VolatileRetentionProfile {
            proportion_const: 75000.0,
            retention_divisor: 40.0,
            greenhouse_candidate: false,
        },
        OrbitalZone::Zone3 => VolatileRetentionProfile {
            proportion_const: 2500.0,
            retention_divisor: 100.0,
            greenhouse_candidate: false,
        },
    }
}

/// Chooses coarse volatile-retention coefficients from the body's broad climate state.
///
/// The atmosphere research pass emphasized that atmospheric outcome is not purely an orbital-zone
/// question. Dry rocky, temperate volatile-bearing, icy, and gas-envelope bodies can occupy very
/// different atmospheric regimes even at similar insolation.
fn volatile_retention_profile_for_body(body: &Body) -> VolatileRetentionProfile {
    match body.climate_class() {
        Some(BodyClimateClass::DryRocky) => VolatileRetentionProfile {
            proportion_const: 100000.0,
            retention_divisor: 10.0,
            greenhouse_candidate: true,
        },
        Some(BodyClimateClass::TemperateRocky) => VolatileRetentionProfile {
            proportion_const: 75000.0,
            retention_divisor: 20.0,
            greenhouse_candidate: true,
        },
        Some(BodyClimateClass::IceRich) => VolatileRetentionProfile {
            proportion_const: 120000.0,
            retention_divisor: 4.0,
            greenhouse_candidate: false,
        },
        Some(BodyClimateClass::GasEnvelope) => VolatileRetentionProfile {
            proportion_const: 160000.0,
            retention_divisor: 2.0,
            greenhouse_candidate: false,
        },
        None => legacy_volatile_retention_profile(body.orbit_zone.clone()),
    }
}

/// Returns whether the body should be treated as greenhouse-prone for inventory scaling.
///
/// This stays intentionally simple: climate-specific profiles opt into greenhouse handling, and
/// the star supplies the greenhouse radius threshold.
fn greenhouse_effect_for_body(body: &Body, star: &Star, profile: VolatileRetentionProfile) -> bool {
    profile.greenhouse_candidate && body.a < star.r_greenhouse
}

/// (Fogg 1985)-style volatile inventory estimate using a body-specific retention profile.
///
/// This keeps the shape of the classic relation while allowing formation/climate-aware tuning.
fn vol_inventory_with_profile(
    mass: f64,
    esc_velocity: f64,
    rms_velocity: f64,
    stellar_mass: f64,
    profile: VolatileRetentionProfile,
    greenhouse_effect: bool,
) -> f64 {
    let velocity_ratio = esc_velocity / rms_velocity;
    if velocity_ratio < env::GAS_RETENTION_THRESHOLD {
        return 0.0;
    }

    let earth_units = mass * consts::SUN_MASS_IN_EARTH_MASSES;
    let temp1 = (profile.proportion_const * earth_units) / stellar_mass;
    let temp2 = about(temp1, 0.2);

    if greenhouse_effect {
        temp2
    } else {
        temp2 / profile.retention_divisor
    }
}

/// Infers coarse water / volatile-ice / primordial-gas fractions for atmosphere synthesis.
///
/// When available, this prefers `body.formation_inputs.condensation_fractions` so atmospheric
/// outcome is tied to disk chemistry and volatile delivery, as suggested by the README notes and
/// the protoplanetary chemistry / planet formation references cited there.
fn inferred_composition_fractions(body: &Body) -> (f64, f64, f64) {
    if let Some(inputs) = body.formation_inputs {
        let fractions = inputs.condensation_fractions;
        return (
            (fractions.water_ice + 0.35 * fractions.volatile_ices).clamp(0.0, 1.0),
            (fractions.volatile_ices + 0.25 * fractions.water_ice).clamp(0.0, 1.0),
            fractions.gas.clamp(0.0, 1.0),
        );
    }

    match body.climate_class() {
        Some(BodyClimateClass::DryRocky) => (0.02, 0.03, 0.01),
        Some(BodyClimateClass::TemperateRocky) => (0.18, 0.07, 0.02),
        Some(BodyClimateClass::IceRich) => (0.45, 0.20, 0.03),
        Some(BodyClimateClass::GasEnvelope) => (0.04, 0.10, 0.60),
        None => match body.orbit_zone {
            OrbitalZone::Zone1 => (
                0.01,
                0.02,
                if body.mass_type == MassType::GasGiant {
                    0.60
                } else {
                    0.02
                },
            ),
            OrbitalZone::Zone2 => (
                0.15,
                0.07,
                if body.mass_type == MassType::GasGiant {
                    0.60
                } else {
                    0.02
                },
            ),
            OrbitalZone::Zone3 => (
                0.45,
                0.20,
                if body.mass_type == MassType::GasGiant {
                    0.60
                } else {
                    0.03
                },
            ),
        },
    }
}

/// Supplies a first-pass albedo used to estimate equilibrium temperature before the richer
/// class-based surface iteration runs.
fn seed_albedo_for_body(body: &Body) -> f64 {
    match body.climate_class() {
        Some(BodyClimateClass::DryRocky) => 0.16,
        Some(BodyClimateClass::TemperateRocky) => 0.28,
        Some(BodyClimateClass::IceRich) => 0.45,
        Some(BodyClimateClass::GasEnvelope) => 0.42,
        None => 0.30,
    }
}

/// Inverse of `pressure`, converting a pressure estimate back into legacy volatile inventory units.
fn inverse_pressure(surf_pressure_mb: f64, equat_radius: f64, grav: f64) -> f64 {
    if surf_pressure_mb <= 0.0 || equat_radius <= 0.0 || grav <= 0.0 {
        return 0.0;
    }

    let equat_radius_ratio = env::KM_EARTH_RADIUS / equat_radius;
    surf_pressure_mb * equat_radius_ratio.powf(2.0) / grav
}

/// Returns a smooth species-by-species retention fraction for a candidate gas.
///
/// A key goal of the atmosphere overhaul was to avoid collapsing the model to a single retained
/// molecule threshold. Evaluating each species separately lets marginal rocky worlds keep some
/// `CO2`, `H2O`, or `Ar` even when they struggle to retain `N2`.
fn species_retention_fraction(body: &Body, star: &Star, radius_km: f64, gas: AtmosphereGas) -> f64 {
    let escape_velocity = escape_vel(body.mass_in_sols, radius_km);
    let rms_velocity = rms_vel(gas.molecular_weight(), body.a, star.luminosity_in_sols);
    if escape_velocity <= 0.0 || rms_velocity <= 0.0 {
        return 0.0;
    }

    // `normalized == 1` corresponds roughly to the classic retention threshold. We deliberately
    // soften the transition around that boundary to avoid brittle all-or-nothing atmospheres.
    let normalized = (escape_velocity / rms_velocity) / env::GAS_RETENTION_THRESHOLD;
    if normalized <= 0.45 {
        0.0
    } else {
        // The sub-linear exponent keeps the transition smooth while still favoring bodies with
        // a comfortable escape-velocity margin over thermal velocity.
        ((normalized - 0.45) / 0.55).clamp(0.0, 1.0).powf(0.85)
    }
}

/// Adds a weighted gas contribution to a composition pool, merging duplicate species.
fn add_gas_weight(weights: &mut Vec<(AtmosphereGas, f64)>, gas: AtmosphereGas, amount: f64) {
    if amount <= 0.0 {
        return;
    }

    if let Some((_, existing)) = weights.iter_mut().find(|(existing_gas, _)| *existing_gas == gas) {
        *existing += amount;
    } else {
        weights.push((gas, amount));
    }
}

/// Distributes a pressure budget across weighted gases, then applies per-species retention.
fn add_retained_pressure(
    contributions: &mut Vec<(AtmosphereGas, f64)>,
    weights: &[(AtmosphereGas, f64)],
    pressure_budget_mb: f64,
    body: &Body,
    star: &Star,
    radius_km: f64,
) {
    if pressure_budget_mb <= 0.0 || weights.is_empty() {
        return;
    }

    let total_weight = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
    if total_weight <= 0.0 {
        return;
    }

    for (gas, weight) in weights {
        let retention = species_retention_fraction(body, star, radius_km, *gas);
        if retention <= 0.0 {
            continue;
        }

        let retained_pressure = pressure_budget_mb * (*weight / total_weight) * retention;
        if retained_pressure <= 0.0 {
            continue;
        }

        if let Some((_, existing)) = contributions.iter_mut().find(|(existing_gas, _)| *existing_gas == *gas) {
            *existing += retained_pressure;
        } else {
            contributions.push((*gas, retained_pressure));
        }
    }
}

/// Looks up the retained fraction of a given species in a synthesized atmosphere.
fn component_fraction(components: &[AtmosphereGasComponent], gas: AtmosphereGas) -> f64 {
    components
        .iter()
        .find(|component| component.gas == gas)
        .map(|component| component.fraction)
        .unwrap_or(0.0)
}

/// Infers broad gameplay/habitability traits from a retained gas mixture.
///
/// The thresholds here are intentionally conservative and readable rather than exhaustive. They
/// are meant to capture broad categories from the atmosphere notes: breathable `N2/O2` bands,
/// greenhouse/steam worlds, reducing atmospheres, and strongly corrosive sulfur/CO2 cases.
fn derive_atmosphere_traits(
    components: &[AtmosphereGasComponent],
    surf_pressure_mb: f64,
    eq_temp_k: f64,
) -> AtmosphereTraits {
    let hydrogen_fraction = component_fraction(components, AtmosphereGas::Hydrogen);
    let helium_fraction = component_fraction(components, AtmosphereGas::Helium);
    let methane_fraction = component_fraction(components, AtmosphereGas::Methane);
    let ammonia_fraction = component_fraction(components, AtmosphereGas::Ammonia);
    let water_vapor_fraction = component_fraction(components, AtmosphereGas::WaterVapor);
    let nitrogen_fraction = component_fraction(components, AtmosphereGas::Nitrogen);
    let oxygen_fraction = component_fraction(components, AtmosphereGas::Oxygen);
    let carbon_dioxide_fraction = component_fraction(components, AtmosphereGas::CarbonDioxide);
    let sulfur_dioxide_fraction = component_fraction(components, AtmosphereGas::SulfurDioxide);
    let inert_fraction = nitrogen_fraction + component_fraction(components, AtmosphereGas::Argon);

    let hydrogen_rich = hydrogen_fraction >= 0.50 || hydrogen_fraction + helium_fraction >= 0.72;
    let steam = water_vapor_fraction >= 0.22 && eq_temp_k >= 330.0 && surf_pressure_mb >= 100.0;
    let oxidizing = oxygen_fraction >= 0.18;
    let reducing = hydrogen_rich || methane_fraction >= 0.05 || ammonia_fraction >= 0.005;
    let corrosive = (sulfur_dioxide_fraction >= 0.02 && surf_pressure_mb >= 150.0)
        || (steam && sulfur_dioxide_fraction >= 0.005)
        || (carbon_dioxide_fraction >= 0.75 && eq_temp_k >= 500.0);
    let toxic = hydrogen_fraction >= 0.70
        || carbon_dioxide_fraction >= 0.10
        || ammonia_fraction >= 0.001
        || sulfur_dioxide_fraction >= 0.001;
    let breathable = (350.0..=3_000.0).contains(&surf_pressure_mb)
        && (0.16..=0.30).contains(&oxygen_fraction)
        && inert_fraction >= 0.45
        && carbon_dioxide_fraction <= 0.03
        && sulfur_dioxide_fraction <= 0.001
        && ammonia_fraction <= 0.0001
        && hydrogen_fraction <= 0.05;
    let tainted = !breathable
        && (toxic
            || steam
            || methane_fraction >= 0.03
            || carbon_dioxide_fraction >= 0.03
            || water_vapor_fraction >= 0.10
            || !(250.0..=320.0).contains(&eq_temp_k)
            || !(250.0..=2_500.0).contains(&surf_pressure_mb));

    AtmosphereTraits {
        breathable,
        tainted,
        toxic,
        oxidizing,
        reducing,
        corrosive,
        steam,
        hydrogen_rich,
    }
}

/// Builds the detailed retained atmosphere mixture for a body.
///
/// Roughly speaking this function:
/// 1. infers bulk volatile / water / gas reservoirs,
/// 2. estimates secondary and primordial pressure budgets,
/// 3. assigns likely gas weights for those reservoirs,
/// 4. applies species-by-species retention,
/// 5. derives qualitative traits and summary metrics.
///
/// The goal is to reflect the broad themes from the README atmosphere notes—volatile delivery,
/// condensation chemistry, secondary outgassing, primordial `H2/He` envelopes, and escape—while
/// still staying lightweight and game-friendly.
fn synthesize_atmosphere_model(
    body: &Body,
    star: &Star,
    radius_km: f64,
    surface_grav_gees: f64,
    profile: VolatileRetentionProfile,
    greenhouse_effect: bool,
) -> AtmosphereModel {
    // Start from formation-aware composition hints so atmosphere generation is not driven solely
    // by present-day orbit or by a single retained-molecule cutoff.
    let (water_inventory, volatile_inventory_fraction, gas_fraction) = inferred_composition_fractions(body);
    let climate_class = body.climate_class();
    let eq_temp_k = eff_temp(
        earth_equivalent_insolation_radius_au(star),
        body.a,
        seed_albedo_for_body(body),
    );
    let earth_mass = body.mass_in_earth_masses().max(0.0);
    let reservoir_multiplier = match climate_class {
        Some(BodyClimateClass::DryRocky) => {
            (0.01 + 0.30 * volatile_inventory_fraction + 0.18 * water_inventory + 0.08 * gas_fraction)
                .clamp(0.005, 0.25)
        }
        Some(BodyClimateClass::TemperateRocky) => {
            (0.10 + 0.85 * volatile_inventory_fraction + 0.65 * water_inventory + 0.25 * gas_fraction).clamp(0.08, 0.95)
        }
        Some(BodyClimateClass::IceRich) => {
            (0.18 + 0.95 * volatile_inventory_fraction + 0.85 * water_inventory + 0.25 * gas_fraction).clamp(0.15, 1.25)
        }
        Some(BodyClimateClass::GasEnvelope) => {
            (0.30 + 0.40 * volatile_inventory_fraction + 0.30 * water_inventory + 2.60 * gas_fraction).clamp(0.40, 3.00)
        }
        None => {
            (0.08 + 0.80 * volatile_inventory_fraction + 0.60 * water_inventory + 1.60 * gas_fraction).clamp(0.05, 1.50)
        }
    };
    // This still uses the classic volatile-inventory shape, but the multiplier lets climate and
    // formation chemistry influence how much source material the atmosphere starts from.
    let raw_volatile_inventory = if earth_mass > 0.0 {
        let base = (profile.proportion_const * earth_mass) / star.mass_in_sols.max(0.1) * reservoir_multiplier;
        if greenhouse_effect {
            base
        } else {
            base / profile.retention_divisor
        }
    } else {
        0.0
    };
    let base_pressure = pressure(raw_volatile_inventory, radius_km, surface_grav_gees);
    // Compress very large raw pressures so volatile-rich worlds still spread out across useful
    // gameplay classes instead of immediately saturating the top end.
    let compressed_base_pressure = if base_pressure > 0.0 {
        env::EARTH_SURF_PRES_IN_MILLIBARS * (base_pressure / env::EARTH_SURF_PRES_IN_MILLIBARS).ln_1p() / 4.0
    } else {
        0.0
    };
    let gravity_factor = 0.85 + 0.15 * surface_grav_gees.clamp(0.1, 4.0);
    let secondary_climate_boost = match climate_class {
        Some(BodyClimateClass::DryRocky) => 0.35,
        Some(BodyClimateClass::TemperateRocky) => 0.95,
        Some(BodyClimateClass::IceRich) => 1.20,
        Some(BodyClimateClass::GasEnvelope) => 0.70,
        None => 0.75,
    };
    let secondary_reservoir_factor = match climate_class {
        Some(BodyClimateClass::DryRocky) => {
            (0.04 + 0.35 * volatile_inventory_fraction + 0.22 * water_inventory).clamp(0.01, 0.35)
        }
        Some(BodyClimateClass::TemperateRocky) => {
            (0.85 + 1.35 * volatile_inventory_fraction + 1.15 * water_inventory).clamp(0.45, 1.60)
        }
        Some(BodyClimateClass::IceRich) => {
            (0.60 + 1.40 * volatile_inventory_fraction + 1.20 * water_inventory).clamp(0.55, 1.90)
        }
        Some(BodyClimateClass::GasEnvelope) => {
            (0.30 + 0.80 * volatile_inventory_fraction + 0.60 * water_inventory).clamp(0.20, 1.20)
        }
        None => (0.25 + 1.00 * volatile_inventory_fraction + 0.85 * water_inventory).clamp(0.10, 1.40),
    };
    let secondary_budget_mb =
        compressed_base_pressure * gravity_factor * secondary_reservoir_factor * secondary_climate_boost;
    let primordial_budget_mb = if matches!(climate_class, Some(BodyClimateClass::GasEnvelope))
        || body.mass_type == MassType::GasGiant
        || gas_fraction >= 0.12
    {
        compressed_base_pressure
            * gravity_factor
            * gas_fraction.clamp(0.02, 1.0)
            * (2.5
                + 8.0 * gas_fraction
                + if matches!(climate_class, Some(BodyClimateClass::GasEnvelope)) {
                    4.0
                } else {
                    0.0
                }
                + if body.mass_type == MassType::GasGiant { 3.0 } else { 0.0 })
    } else {
        0.0
    };
    let hot_volatile = (water_inventory >= 0.10 || volatile_inventory_fraction >= 0.15) && eq_temp_k >= 360.0;
    // Hot volatile-bearing cases get a pressure floor so steam/runaway worlds do not end up
    // looking implausibly thin after the compact class mapping.
    let thermal_pressure_floor = if hot_volatile {
        let reservoir_strength = (0.7 * water_inventory + volatile_inventory_fraction).clamp(0.0, 1.25);
        let heat_term = ((eq_temp_k - 330.0) / 220.0).clamp(0.0, 2.0);
        env::EARTH_SURF_PRES_IN_MILLIBARS
            * reservoir_strength
            * (1.2 + 1.4 * heat_term + if greenhouse_effect { 0.8 } else { 0.0 })
    } else {
        0.0
    };

    // Secondary weights represent atmospheres supplied by outgassing, volatile release, and
    // surface-atmosphere cycling.
    let mut secondary_weights = Vec::new();
    match climate_class {
        Some(BodyClimateClass::DryRocky) => {
            add_gas_weight(&mut secondary_weights, AtmosphereGas::CarbonDioxide, 0.50);
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Nitrogen, 0.18);
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Argon, 0.04);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::WaterVapor,
                0.02 + 0.18 * water_inventory,
            );
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::SulfurDioxide,
                0.005
                    + if eq_temp_k >= 450.0 {
                        0.08 * volatile_inventory_fraction
                    } else {
                        0.01 * volatile_inventory_fraction
                    },
            );
        }
        Some(BodyClimateClass::TemperateRocky) => {
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Nitrogen, 0.46);
            add_gas_weight(&mut secondary_weights, AtmosphereGas::CarbonDioxide, 0.16);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::WaterVapor,
                0.08 + 0.15 * water_inventory,
            );
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Argon, 0.03);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::Oxygen,
                0.06 + 0.10 * (water_inventory * (1.0 - ((eq_temp_k - 290.0).abs() / 140.0)).clamp(0.0, 1.0)),
            );
            if eq_temp_k < 245.0 {
                add_gas_weight(
                    &mut secondary_weights,
                    AtmosphereGas::Methane,
                    0.03 * volatile_inventory_fraction,
                );
            }
        }
        Some(BodyClimateClass::IceRich) => {
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Nitrogen, 0.24);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::WaterVapor,
                0.12 + 0.25 * water_inventory,
            );
            add_gas_weight(&mut secondary_weights, AtmosphereGas::CarbonDioxide, 0.18);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::Methane,
                0.12 * volatile_inventory_fraction.clamp(0.0, 1.0),
            );
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::Ammonia,
                0.07 * volatile_inventory_fraction,
            );
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Argon, 0.03);
        }
        Some(BodyClimateClass::GasEnvelope) => {
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::WaterVapor,
                0.04 + 0.06 * water_inventory,
            );
            if eq_temp_k < 280.0 {
                add_gas_weight(
                    &mut secondary_weights,
                    AtmosphereGas::Methane,
                    0.02 + 0.05 * volatile_inventory_fraction,
                );
            }
            if eq_temp_k < 220.0 {
                add_gas_weight(
                    &mut secondary_weights,
                    AtmosphereGas::Ammonia,
                    0.01 + 0.03 * volatile_inventory_fraction,
                );
            }
            if eq_temp_k >= 250.0 {
                add_gas_weight(&mut secondary_weights, AtmosphereGas::CarbonDioxide, 0.02);
            }
        }
        None => {
            add_gas_weight(&mut secondary_weights, AtmosphereGas::Nitrogen, 0.20);
            add_gas_weight(&mut secondary_weights, AtmosphereGas::CarbonDioxide, 0.20);
            add_gas_weight(
                &mut secondary_weights,
                AtmosphereGas::WaterVapor,
                0.05 + 0.10 * water_inventory,
            );
        }
    }

    add_gas_weight(
        &mut secondary_weights,
        AtmosphereGas::WaterVapor,
        water_inventory * if hot_volatile { 0.45 } else { 0.20 },
    );
    add_gas_weight(
        &mut secondary_weights,
        AtmosphereGas::CarbonDioxide,
        volatile_inventory_fraction * if eq_temp_k >= 380.0 { 0.30 } else { 0.18 },
    );
    add_gas_weight(
        &mut secondary_weights,
        AtmosphereGas::Nitrogen,
        volatile_inventory_fraction * 0.12,
    );
    add_gas_weight(
        &mut secondary_weights,
        AtmosphereGas::Argon,
        volatile_inventory_fraction * 0.02,
    );
    if eq_temp_k >= 500.0 {
        add_gas_weight(
            &mut secondary_weights,
            AtmosphereGas::SulfurDioxide,
            volatile_inventory_fraction * 0.10,
        );
    } else if eq_temp_k >= 380.0 {
        add_gas_weight(
            &mut secondary_weights,
            AtmosphereGas::SulfurDioxide,
            volatile_inventory_fraction * 0.04,
        );
    }
    if water_inventory >= 0.08 && eq_temp_k >= 340.0 {
        add_gas_weight(
            &mut secondary_weights,
            AtmosphereGas::Oxygen,
            0.12 * ((eq_temp_k - 340.0) / 220.0).clamp(0.0, 1.0),
        );
    }

    // Primordial weights represent retained nebular gas envelopes on gas-rich worlds.
    let mut primordial_weights = Vec::new();
    if primordial_budget_mb > 0.0 {
        add_gas_weight(
            &mut primordial_weights,
            AtmosphereGas::Hydrogen,
            0.72 + 0.40 * gas_fraction,
        );
        add_gas_weight(
            &mut primordial_weights,
            AtmosphereGas::Helium,
            0.18 + 0.10 * gas_fraction,
        );
        if eq_temp_k < 280.0 {
            add_gas_weight(
                &mut primordial_weights,
                AtmosphereGas::Methane,
                0.02 + 0.04 * volatile_inventory_fraction,
            );
        }
        if eq_temp_k < 220.0 {
            add_gas_weight(
                &mut primordial_weights,
                AtmosphereGas::Ammonia,
                0.01 + 0.03 * volatile_inventory_fraction,
            );
        }
        add_gas_weight(
            &mut primordial_weights,
            AtmosphereGas::WaterVapor,
            if eq_temp_k > 280.0 {
                0.04 * water_inventory
            } else {
                0.01 * water_inventory
            },
        );
        if eq_temp_k >= 300.0 {
            add_gas_weight(&mut primordial_weights, AtmosphereGas::CarbonDioxide, 0.01);
        }
    }

    // Convert reservoir budgets into a retained atmosphere by applying species-specific survival.
    let mut retained_pressure_by_gas = Vec::new();
    add_retained_pressure(
        &mut retained_pressure_by_gas,
        &secondary_weights,
        secondary_budget_mb,
        body,
        star,
        radius_km,
    );
    add_retained_pressure(
        &mut retained_pressure_by_gas,
        &primordial_weights,
        primordial_budget_mb,
        body,
        star,
        radius_km,
    );

    let mut total_pressure_mb = retained_pressure_by_gas
        .iter()
        .map(|(_, pressure_mb)| *pressure_mb)
        .sum::<f64>();
    if thermal_pressure_floor > total_pressure_mb {
        let delta = thermal_pressure_floor - total_pressure_mb;
        // Inject the floor as a hot volatile blend rather than a single gas so downstream
        // classification and UI still see a meaningful composition.
        add_gas_weight(&mut retained_pressure_by_gas, AtmosphereGas::WaterVapor, delta * 0.65);
        add_gas_weight(
            &mut retained_pressure_by_gas,
            AtmosphereGas::CarbonDioxide,
            delta * if eq_temp_k >= 450.0 { 0.25 } else { 0.20 },
        );
        if eq_temp_k >= 500.0 {
            add_gas_weight(
                &mut retained_pressure_by_gas,
                AtmosphereGas::SulfurDioxide,
                delta * 0.10,
            );
        } else {
            add_gas_weight(&mut retained_pressure_by_gas, AtmosphereGas::Nitrogen, delta * 0.10);
            add_gas_weight(&mut retained_pressure_by_gas, AtmosphereGas::Oxygen, delta * 0.05);
        }
        total_pressure_mb = retained_pressure_by_gas
            .iter()
            .map(|(_, pressure_mb)| *pressure_mb)
            .sum::<f64>();
    }

    // If nothing survives retention, return an effectively airless model while still preserving
    // formation-derived reservoir hints for downstream temperature/surface logic.
    if total_pressure_mb <= 0.0 {
        return AtmosphereModel {
            water_inventory,
            volatile_inventory_fraction,
            gas_fraction,
            eq_temp_k,
            smallest_mw_retained: molecule_limit(body.mass_in_sols, radius_km),
            ..Default::default()
        };
    }

    retained_pressure_by_gas.sort_by(|a, b| b.1.total_cmp(&a.1));

    // Convert retained pressures into stable components and summary statistics used by the UI,
    // classification layer, and tests.
    let volatile_gas_inventory = inverse_pressure(total_pressure_mb, radius_km, surface_grav_gees);
    let mut components = retained_pressure_by_gas
        .iter()
        .map(|(gas, pressure_mb)| AtmosphereGasComponent {
            gas: *gas,
            fraction: (*pressure_mb / total_pressure_mb).clamp(0.0, 1.0),
            partial_pressure_mb: *pressure_mb,
            retained_inventory: volatile_gas_inventory * (*pressure_mb / total_pressure_mb).clamp(0.0, 1.0),
        })
        .collect::<Vec<_>>();
    components.retain(|component| component.partial_pressure_mb > 0.01);

    let dominant_gas = components.first().map(|component| component.gas);
    let mean_molecular_weight = components
        .iter()
        .map(|component| component.gas.molecular_weight() * component.fraction)
        .sum::<f64>();
    let smallest_mw_retained = components
        .iter()
        .filter(|component| component.fraction >= 0.02 || component.partial_pressure_mb >= 10.0)
        .map(|component| component.gas.molecular_weight())
        .min_by(|a, b| a.total_cmp(b))
        .or_else(|| dominant_gas.map(|gas| gas.molecular_weight()))
        .unwrap_or_else(|| molecule_limit(body.mass_in_sols, radius_km));
    let traits = derive_atmosphere_traits(&components, total_pressure_mb, eq_temp_k);

    AtmosphereModel {
        components,
        dominant_gas,
        traits,
        volatile_gas_inventory,
        surf_pressure_mb: total_pressure_mb,
        smallest_mw_retained,
        mean_molecular_weight,
        water_inventory,
        volatile_inventory_fraction,
        gas_fraction,
        eq_temp_k,
    }
}

/// Projects the detailed atmosphere model into the smaller feature set used by the classifier.
fn build_atmosphere_signals(
    body: &Body,
    atmosphere_model: &AtmosphereModel,
    surface_grav_gees: f64,
    greenhouse_effect: bool,
) -> AtmosphereSignals {
    let climate_class = body.climate_class();

    AtmosphereSignals {
        eq_temp_k: atmosphere_model.eq_temp_k,
        smallest_mw_retained: atmosphere_model.smallest_mw_retained,
        mean_molecular_weight: atmosphere_model.mean_molecular_weight,
        volatile_gas_inventory: atmosphere_model.volatile_gas_inventory,
        surface_grav_gees,
        pressure_proxy_mb: atmosphere_model.surf_pressure_mb,
        mass_earth: body.mass_in_earth_masses(),
        gas_fraction: atmosphere_model.gas_fraction,
        water_inventory: atmosphere_model.water_inventory,
        volatile_inventory_fraction: atmosphere_model.volatile_inventory_fraction,
        climate_class,
        greenhouse_effect,
        mass_type: body.mass_type,
        dominant_gas: atmosphere_model.dominant_gas,
        hydrogen_fraction: component_fraction(&atmosphere_model.components, AtmosphereGas::Hydrogen),
        water_vapor_fraction: component_fraction(&atmosphere_model.components, AtmosphereGas::WaterVapor),
        oxygen_fraction: component_fraction(&atmosphere_model.components, AtmosphereGas::Oxygen),
        carbon_dioxide_fraction: component_fraction(&atmosphere_model.components, AtmosphereGas::CarbonDioxide),
        sulfur_dioxide_fraction: component_fraction(&atmosphere_model.components, AtmosphereGas::SulfurDioxide),
        traits: atmosphere_model.traits,
    }
}

/// Maps the retained gas mixture back to the compact `AtmosphereClass` presentation enum.
///
/// Ordering matters here: special cases like `GasGiant`, `H2Rich`, `Steam`, and `Corrosive`
/// should win before the generic pressure-band classes do.
fn classify_atmosphere(signals: AtmosphereSignals) -> AtmosphereClass {
    let pressure_atm_proxy = signals.pressure_proxy_mb / env::EARTH_SURF_PRES_IN_MILLIBARS;
    let breathable_mw = (18.0..=42.0).contains(&signals.mean_molecular_weight);
    let temperate_world = matches!(signals.climate_class, Some(BodyClimateClass::TemperateRocky));
    let breathable_band = signals.traits.breathable || breathable_mw || temperate_world;
    let volatile_rich = signals.volatile_inventory_fraction >= 0.12 || signals.volatile_gas_inventory >= 700.0;
    let hot_volatile = signals.traits.steam
        || ((signals.water_inventory >= 0.10 || signals.volatile_inventory_fraction >= 0.15)
            && signals.eq_temp_k >= 360.0);
    let gas_rich = matches!(signals.climate_class, Some(BodyClimateClass::GasEnvelope))
        || signals.gas_fraction >= 0.18
        || signals.hydrogen_fraction >= 0.45;
    let tainted_candidate = signals.traits.tainted
        || signals.traits.toxic
        || matches!(signals.climate_class, Some(BodyClimateClass::DryRocky))
        || signals.eq_temp_k > 325.0
        || signals.eq_temp_k < 235.0
        || (signals.oxygen_fraction > 0.35 && pressure_atm_proxy >= 0.5)
        || (!breathable_mw && !temperate_world)
        || signals.greenhouse_effect;

    // Massive or strongly hydrogen-dominated bodies stay in the gas-giant family even if their
    // pressure proxy could also fit inside one of the lower bands.
    if signals.mass_type == MassType::GasGiant
        || signals.mass_earth >= 30.0
        || (signals.hydrogen_fraction >= 0.85 && pressure_atm_proxy >= 8.0)
        || (gas_rich && pressure_atm_proxy >= 12.0)
    {
        return AtmosphereClass::GasGiant;
    }

    if (signals.surface_grav_gees > 2.0 && pressure_atm_proxy < 0.12)
        || (signals.surface_grav_gees < 0.25 && pressure_atm_proxy > 4.0)
    {
        return AtmosphereClass::Unusual;
    }

    if signals.traits.hydrogen_rich && signals.mass_earth >= 2.0 && pressure_atm_proxy >= 0.8 {
        return AtmosphereClass::H2Rich;
    }

    if signals.traits.corrosive && signals.eq_temp_k >= 850.0 && (pressure_atm_proxy >= 2.0 || hot_volatile) {
        return AtmosphereClass::Insidious;
    }

    if signals.traits.corrosive && (pressure_atm_proxy >= 0.15 || hot_volatile)
        || (signals.eq_temp_k >= 500.0
            && (pressure_atm_proxy >= 1.0
                || (signals.carbon_dioxide_fraction >= 0.75 && pressure_atm_proxy >= 0.5)
                || signals.sulfur_dioxide_fraction >= 0.02
                || (hot_volatile && signals.greenhouse_effect)))
    {
        return AtmosphereClass::Corrosive;
    }

    if signals.traits.steam && (volatile_rich || signals.water_vapor_fraction >= 0.22) && pressure_atm_proxy >= 0.3 {
        return AtmosphereClass::Steam;
    }

    if signals.pressure_proxy_mb < 1.0 || signals.volatile_gas_inventory <= 0.01 {
        return AtmosphereClass::Airless;
    }

    if pressure_atm_proxy < 0.09 {
        return AtmosphereClass::Trace;
    }

    if signals.surface_grav_gees >= 1.2 && breathable_band && pressure_atm_proxy <= 0.50 {
        return AtmosphereClass::ThinLow;
    }

    if signals.surface_grav_gees >= 1.2 && breathable_band && pressure_atm_proxy >= 2.50 {
        return AtmosphereClass::DenseHigh;
    }

    if gas_rich && pressure_atm_proxy >= 0.10 {
        return AtmosphereClass::Exotic;
    }

    if !breathable_band
        && pressure_atm_proxy >= 0.10
        && (signals.smallest_mw_retained < 4.0
            || signals.hydrogen_fraction >= 0.15
            || signals.mean_molecular_weight < 8.0
            || signals.mean_molecular_weight > 60.0
            || signals.eq_temp_k < 220.0
            || signals.dominant_gas == Some(AtmosphereGas::Methane)
            || signals.dominant_gas == Some(AtmosphereGas::Ammonia))
    {
        return AtmosphereClass::Exotic;
    }

    if pressure_atm_proxy < 0.42 {
        return if tainted_candidate {
            AtmosphereClass::VeryThinTainted
        } else {
            AtmosphereClass::VeryThin
        };
    }

    if pressure_atm_proxy < 0.70 {
        return if tainted_candidate {
            AtmosphereClass::ThinTainted
        } else {
            AtmosphereClass::Thin
        };
    }

    if pressure_atm_proxy < 1.49 {
        return if tainted_candidate {
            AtmosphereClass::StandardTainted
        } else {
            AtmosphereClass::Standard
        };
    }

    if pressure_atm_proxy < 2.49 {
        return if tainted_candidate {
            AtmosphereClass::DenseTainted
        } else {
            AtmosphereClass::Dense
        };
    }

    AtmosphereClass::DenseHigh
}

/// Canonical pressure bands for the compact atmosphere presentation classes.
fn atmosphere_pressure_range_mb(class: AtmosphereClass) -> (f64, f64) {
    match class {
        AtmosphereClass::Airless => (0.0, 0.0),
        AtmosphereClass::Trace => (1.0, 90.0),
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => (100.0, 420.0),
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted => (430.0, 700.0),
        AtmosphereClass::Standard | AtmosphereClass::StandardTainted => (710.0, 1490.0),
        AtmosphereClass::Dense | AtmosphereClass::DenseTainted => (1500.0, 2490.0),
        AtmosphereClass::Exotic => (200.0, 4000.0),
        AtmosphereClass::Corrosive => (1500.0, 12000.0),
        AtmosphereClass::Insidious => (3000.0, 30000.0),
        AtmosphereClass::DenseHigh => (2500.0, 6000.0),
        AtmosphereClass::ThinLow => (150.0, 500.0),
        AtmosphereClass::Steam => (2000.0, 20000.0),
        AtmosphereClass::H2Rich => (3000.0, 50000.0),
        AtmosphereClass::GasGiant => (15000.0, 250000.0),
        AtmosphereClass::Unusual => (50.0, 15000.0),
    }
}

/// Clamps the synthesized pressure proxy into the canonical band for the chosen class.
fn class_surface_pressure_mb(class: AtmosphereClass, signals: AtmosphereSignals) -> f64 {
    if matches!(class, AtmosphereClass::Airless) {
        return 0.0;
    }

    let (min_mb, max_mb) = atmosphere_pressure_range_mb(class);
    let adjusted_pressure = signals.pressure_proxy_mb
        * match class {
            AtmosphereClass::Trace => 1.0,
            AtmosphereClass::Steam => 1.05,
            AtmosphereClass::Corrosive => 1.08,
            AtmosphereClass::Insidious => 1.12,
            AtmosphereClass::H2Rich => 1.10,
            AtmosphereClass::GasGiant => 1.15,
            _ => 1.0,
        };

    adjusted_pressure.clamp(min_mb, max_mb)
}

/// First-pass albedo guess associated with the compact atmosphere class.
fn base_albedo_for_class(class: AtmosphereClass, climate_class: Option<BodyClimateClass>) -> f64 {
    let mut base: f64 = match class {
        AtmosphereClass::Airless => 0.12,
        AtmosphereClass::Trace => 0.15,
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => 0.18,
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted | AtmosphereClass::ThinLow => 0.20,
        AtmosphereClass::Standard | AtmosphereClass::StandardTainted => 0.24,
        AtmosphereClass::Dense | AtmosphereClass::DenseTainted => 0.28,
        AtmosphereClass::Exotic => 0.26,
        AtmosphereClass::Corrosive => 0.42,
        AtmosphereClass::Insidious => 0.50,
        AtmosphereClass::DenseHigh => 0.32,
        AtmosphereClass::Steam => 0.55,
        AtmosphereClass::H2Rich => 0.45,
        AtmosphereClass::GasGiant => 0.50,
        AtmosphereClass::Unusual => 0.30,
    };

    if matches!(climate_class, Some(BodyClimateClass::IceRich))
        && matches!(
            class,
            AtmosphereClass::Airless
                | AtmosphereClass::Trace
                | AtmosphereClass::VeryThin
                | AtmosphereClass::VeryThinTainted
        )
    {
        base += 0.12;
    }

    base.clamp(0.05, 0.75)
}

/// Returns an extra greenhouse warming term associated with the classified atmosphere.
///
/// This is intentionally heuristic: the detailed gas mix determines class/traits first, then this
/// helper turns that class into a stable warming term for the surface iteration.
fn greenhouse_adjustment_k(
    class: AtmosphereClass,
    surf_pressure_mb: f64,
    eq_temp_k: f64,
    water_inventory: f64,
    greenhouse_effect: bool,
) -> f64 {
    if matches!(class, AtmosphereClass::Airless) || eq_temp_k <= 0.0 {
        return 0.0;
    }

    let pressure_term = (surf_pressure_mb / env::EARTH_SURF_PRES_IN_MILLIBARS).max(0.0).ln_1p();
    let water_term = water_inventory.clamp(0.0, 1.0);
    let greenhouse_bonus = if greenhouse_effect { 18.0 } else { 0.0 };
    let base = match class {
        AtmosphereClass::Airless => 0.0,
        AtmosphereClass::Trace => 2.0 + 3.0 * pressure_term,
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => 5.0 + 6.0 * pressure_term,
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted | AtmosphereClass::ThinLow => {
            9.0 + 10.0 * pressure_term + 4.0 * water_term
        }
        AtmosphereClass::Standard | AtmosphereClass::StandardTainted => 14.0 + 12.0 * pressure_term + 8.0 * water_term,
        AtmosphereClass::Dense | AtmosphereClass::DenseTainted => 20.0 + 16.0 * pressure_term + 10.0 * water_term,
        AtmosphereClass::Exotic => 18.0 + 14.0 * pressure_term,
        AtmosphereClass::Corrosive => 60.0 + 24.0 * pressure_term + greenhouse_bonus,
        AtmosphereClass::Insidious => 110.0 + 30.0 * pressure_term + greenhouse_bonus,
        AtmosphereClass::DenseHigh => 28.0 + 18.0 * pressure_term + 8.0 * water_term,
        AtmosphereClass::Steam => 80.0 + 20.0 * pressure_term + 18.0 * water_term + greenhouse_bonus,
        AtmosphereClass::H2Rich => 65.0 + 16.0 * pressure_term,
        AtmosphereClass::GasGiant => 120.0 + 24.0 * pressure_term,
        AtmosphereClass::Unusual => 22.0 + 18.0 * pressure_term,
    };

    if eq_temp_k < 120.0 {
        base * 0.4
    } else {
        base
    }
}

/// Estimates how much of the body's water inventory can exist at the surface.
fn surface_water_fraction(
    class: AtmosphereClass,
    climate_class: Option<BodyClimateClass>,
    water_inventory: f64,
    surf_temp_k: f64,
    boil_point_k: f64,
    surf_pressure_mb: f64,
) -> f64 {
    let climate_bias = match climate_class {
        Some(BodyClimateClass::DryRocky) => 0.20,
        Some(BodyClimateClass::TemperateRocky) => 0.85,
        Some(BodyClimateClass::IceRich) => 1.15,
        Some(BodyClimateClass::GasEnvelope) => 0.0,
        None => 0.60,
    };
    let class_bias = match class {
        AtmosphereClass::Airless => 0.35,
        AtmosphereClass::Trace => 0.45,
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => 0.55,
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted | AtmosphereClass::ThinLow => 0.72,
        AtmosphereClass::Standard | AtmosphereClass::StandardTainted => 0.95,
        AtmosphereClass::Dense | AtmosphereClass::DenseTainted => 0.90,
        AtmosphereClass::DenseHigh => 0.75,
        AtmosphereClass::Exotic => 0.55,
        AtmosphereClass::Corrosive => 0.08,
        AtmosphereClass::Insidious => 0.03,
        AtmosphereClass::Steam => 0.12,
        AtmosphereClass::H2Rich | AtmosphereClass::GasGiant => 0.0,
        AtmosphereClass::Unusual => 0.60,
    };

    let mut surface_water = (water_inventory * climate_bias * class_bias).clamp(0.0, 1.0);

    if surface_water == 0.0 {
        return 0.0;
    }

    if surf_pressure_mb < 6.1 && surf_temp_k > env::FREEZING_POINT_OF_WATER + 15.0 {
        surface_water *= 0.15;
    }

    if boil_point_k > 0.0 && surf_temp_k >= boil_point_k - 5.0 {
        surface_water *= if matches!(class, AtmosphereClass::Steam) {
            0.20
        } else {
            0.05
        };
    }

    if surf_temp_k < 170.0 {
        surface_water = surface_water.max((water_inventory * 0.50).min(1.0));
    }

    surface_water.clamp(0.0, 1.0)
}

/// Returns plausible cloud-cover bounds for a given compact atmosphere class.
fn cloud_cover_range(class: AtmosphereClass) -> (f64, f64) {
    match class {
        AtmosphereClass::Airless => (0.0, 0.0),
        AtmosphereClass::Trace => (0.0, 0.08),
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => (0.02, 0.15),
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted | AtmosphereClass::ThinLow => (0.05, 0.40),
        AtmosphereClass::Standard | AtmosphereClass::StandardTainted => (0.10, 0.75),
        AtmosphereClass::Dense | AtmosphereClass::DenseTainted => (0.25, 0.90),
        AtmosphereClass::Exotic => (0.15, 0.80),
        AtmosphereClass::Corrosive => (0.50, 1.00),
        AtmosphereClass::Insidious => (0.65, 1.00),
        AtmosphereClass::DenseHigh => (0.30, 0.95),
        AtmosphereClass::Steam => (0.70, 1.00),
        AtmosphereClass::H2Rich => (0.55, 0.95),
        AtmosphereClass::GasGiant => (0.75, 1.00),
        AtmosphereClass::Unusual => (0.20, 0.90),
    }
}

/// Estimates cloud cover from atmosphere class, pressure, temperature, and water supply.
fn class_cloud_cover(
    class: AtmosphereClass,
    hydrosphere: f64,
    water_inventory: f64,
    surf_temp_k: f64,
    surf_pressure_mb: f64,
) -> f64 {
    let (floor, ceiling) = cloud_cover_range(class);
    if ceiling == 0.0 {
        return 0.0;
    }

    let water_term = hydrosphere.max(water_inventory * 0.35).clamp(0.0, 1.0);
    let humidity_term = (1.0 - ((surf_temp_k - 290.0).abs() / 180.0)).clamp(0.0, 1.0);
    let pressure_term = ((surf_pressure_mb / env::EARTH_SURF_PRES_IN_MILLIBARS).max(0.0).ln_1p() / 3.0).clamp(0.0, 1.0);
    let driver = match class {
        AtmosphereClass::Steam => 0.85,
        AtmosphereClass::H2Rich | AtmosphereClass::GasGiant => 0.70,
        AtmosphereClass::Corrosive | AtmosphereClass::Insidious => 0.78,
        _ => (0.55 * water_term + 0.25 * humidity_term + 0.20 * pressure_term).clamp(0.0, 1.0),
    };

    if water_term <= 0.01
        && !matches!(
            class,
            AtmosphereClass::H2Rich
                | AtmosphereClass::GasGiant
                | AtmosphereClass::Steam
                | AtmosphereClass::Corrosive
                | AtmosphereClass::Insidious
                | AtmosphereClass::Exotic
        )
    {
        return floor.min(0.05);
    }

    (floor + (ceiling - floor) * driver).clamp(0.0, 1.0)
}

/// Estimates how much of the surface water inventory is frozen under the classified climate.
fn class_ice_cover(
    class: AtmosphereClass,
    climate_class: Option<BodyClimateClass>,
    hydrosphere: f64,
    surf_temp_k: f64,
    surf_pressure_mb: f64,
) -> f64 {
    if hydrosphere <= 0.0 {
        return 0.0;
    }

    if matches!(
        class,
        AtmosphereClass::Steam
            | AtmosphereClass::Corrosive
            | AtmosphereClass::Insidious
            | AtmosphereClass::H2Rich
            | AtmosphereClass::GasGiant
    ) {
        return 0.0;
    }

    if surf_temp_k < 160.0 {
        return hydrosphere;
    }

    if surf_pressure_mb < 6.1 {
        return if surf_temp_k < 260.0 {
            hydrosphere
        } else {
            hydrosphere * 0.30
        };
    }

    let climate_bias = if matches!(climate_class, Some(BodyClimateClass::IceRich)) {
        0.35
    } else {
        0.0
    };
    let cold_factor = ((283.0 - surf_temp_k) / 80.0).clamp(0.0, 1.0).powf(1.5);
    (hydrosphere * (cold_factor + climate_bias).clamp(0.0, 1.0)).clamp(0.0, hydrosphere)
}

/// Combines class base albedo with water, cloud, and ice state into the final albedo.
fn class_albedo(
    class: AtmosphereClass,
    climate_class: Option<BodyClimateClass>,
    hydrosphere: f64,
    cloud_cover: f64,
    ice_cover: f64,
) -> f64 {
    let liquid_water = (hydrosphere - ice_cover).max(0.0);
    (base_albedo_for_class(class, climate_class) + 0.20 * cloud_cover + 0.18 * ice_cover + 0.03 * hydrosphere
        - 0.04 * liquid_water)
        .clamp(0.05, 0.85)
}

/// Applies the compact atmosphere class to the legacy surface-temperature and hydrology model.
///
/// This is a short fixed-point iteration: albedo affects temperature, temperature affects
/// water/cloud/ice cover, and those feeds back into albedo.
fn apply_classified_surface_model(
    props: &mut EnviroProperties,
    atmosphere_class: AtmosphereClass,
    climate_class: Option<BodyClimateClass>,
    water_inventory: f64,
    greenhouse_effect: bool,
    earth_equivalent_insolation_radius_au: f64,
) {
    let mut albedo = base_albedo_for_class(atmosphere_class, climate_class);
    let mut hydrosphere = 0.0;
    let mut cloud_cover = 0.0;
    let mut ice_cover = 0.0;
    let mut effective_temp = eff_temp(earth_equivalent_insolation_radius_au, props.a, albedo);
    let mut greenhouse_rise = 0.0;
    let mut surf_temp = effective_temp;

    // Two passes are enough to stabilize the broad feedback loop without turning this helper into
    // a heavyweight climate solver.
    for _ in 0..2 {
        effective_temp = eff_temp(earth_equivalent_insolation_radius_au, props.a, albedo);
        greenhouse_rise = greenhouse_adjustment_k(
            atmosphere_class,
            props.surf_pressure,
            effective_temp,
            water_inventory,
            greenhouse_effect,
        );
        surf_temp = effective_temp + greenhouse_rise;
        hydrosphere = surface_water_fraction(
            atmosphere_class,
            climate_class,
            water_inventory,
            surf_temp,
            props.boil_point,
            props.surf_pressure,
        );
        ice_cover = class_ice_cover(
            atmosphere_class,
            climate_class,
            hydrosphere,
            surf_temp,
            props.surf_pressure,
        );
        cloud_cover = class_cloud_cover(
            atmosphere_class,
            hydrosphere,
            water_inventory,
            surf_temp,
            props.surf_pressure,
        );
        albedo = class_albedo(atmosphere_class, climate_class, hydrosphere, cloud_cover, ice_cover);
    }

    props.hydrosphere = hydrosphere;
    props.cloud_cover = cloud_cover;
    props.ice_cover = ice_cover.min(hydrosphere);
    props.albedo = albedo;
    props.effective_temp_k = effective_temp.max(0.0);
    props.greenhouse_rise_k = greenhouse_rise.max(0.0);
    props.surf_temp = surf_temp.max(0.0);
    props.atmosphere_class = atmosphere_class;
}

/// Computes a body's environmental properties from its orbital/physical parameters and its star.
///
/// For the atmosphere path this layers a modernized synthesis stage on top of the classic
/// `enviro` helpers:
/// 1. infer volatile-retention and greenhouse context,
/// 2. synthesize a retained gas mixture,
/// 3. classify that mixture into a compact `AtmosphereClass`,
/// 4. run the surface/albedo/water iteration using that class.
///
/// See `starform-rust/README.md` (`Atmosphere modeling notes`) for the rationale and paper list.
pub fn compute_enviro_properties_for_body(body: &Body, star: &Star) -> EnviroProperties {
    let climate_profile = volatile_retention_profile_for_body(body);

    // Some bodies may not have density/radius populated (or may be placeholders).
    // Use best-effort fallbacks so the UI can still show something.
    let density_gcc = if body.density_in_grams_per_cc > 0.0 {
        body.density_in_grams_per_cc
    } else {
        empirical_density(body.mass_in_sols, body.a, body.mass_type, star.luminosity_in_sols)
    };

    let radius_km = if body.radius_in_km > 0.0 {
        body.radius_in_km
    } else {
        volume_radius(body.mass_in_sols, density_gcc)
    };

    let smallest_mw_retained = molecule_limit(body.mass_in_sols, radius_km);
    let greenhouse_effect = greenhouse_effect_for_body(body, star, climate_profile);

    let surface_grav_gees = gravity(accel(body.mass_in_sols, radius_km));
    let atmosphere_model = synthesize_atmosphere_model(
        body,
        star,
        radius_km,
        surface_grav_gees,
        climate_profile,
        greenhouse_effect,
    );
    let signals = build_atmosphere_signals(body, &atmosphere_model, surface_grav_gees, greenhouse_effect);
    let atmosphere_class = classify_atmosphere(signals);
    let surf_pressure = class_surface_pressure_mb(atmosphere_class, signals);
    let boil_point = boiling_point(surf_pressure);
    let volatile_gas_inventory = inverse_pressure(surf_pressure, radius_km, surface_grav_gees);
    let atmosphere_components = if matches!(atmosphere_class, AtmosphereClass::Airless) || surf_pressure <= 0.0 {
        Vec::new()
    } else {
        atmosphere_model
            .components
            .iter()
            .map(|component| AtmosphereGasComponent {
                gas: component.gas,
                fraction: component.fraction,
                partial_pressure_mb: component.fraction * surf_pressure,
                retained_inventory: component.fraction * volatile_gas_inventory,
            })
            .collect::<Vec<_>>()
    };
    let dominant_gas = atmosphere_components.first().map(|component| component.gas);
    let atmosphere_traits = if atmosphere_components.is_empty() {
        AtmosphereTraits::default()
    } else {
        derive_atmosphere_traits(&atmosphere_components, surf_pressure, atmosphere_model.eq_temp_k)
    };

    let mut props = EnviroProperties {
        a: body.a,
        radius_in_km: radius_km,
        molec_weight: if atmosphere_components.is_empty() {
            smallest_mw_retained
        } else {
            atmosphere_model.smallest_mw_retained
        },
        mean_molecular_weight: if atmosphere_components.is_empty() {
            0.0
        } else {
            atmosphere_model.mean_molecular_weight
        },
        surf_pressure,
        volatile_gas_inventory,
        boil_point,
        atmosphere_class,
        atmosphere_components,
        dominant_gas,
        atmosphere_traits,
        ..Default::default()
    };

    apply_classified_surface_model(
        &mut props,
        atmosphere_class,
        signals.climate_class,
        signals.water_inventory,
        greenhouse_effect,
        earth_equivalent_insolation_radius_au(star),
    );
    props
}

/// This function, given the orbital radius of a planet in AU, returns the
/// orbital `zone` of the particle.
pub fn orb_zone(orb_radius: f64, stell_luminosity_ratio: f64) -> OrbitalZone {
    if orb_radius < (4.0 * stell_luminosity_ratio.sqrt()) {
        OrbitalZone::Zone1
    } else if orb_radius < (15.0 * stell_luminosity_ratio.sqrt()) {
        OrbitalZone::Zone2
    } else {
        OrbitalZone::Zone3
    }
}

/// The mass is in units of solar masses, and the density is in units of grams/cc.
/// The radius returned is in units of km.
pub fn volume_radius(mass: f64, density: f64) -> f64 {
    let mass_in_grams = mass * consts::SOLAR_MASS_IN_GRAMS;
    let volume = mass_in_grams / density;
    ((3.0 * volume) / (4.0 * std::f64::consts::PI)).powf(1.0 / 3.0) / consts::CM_PER_KM
}

/// Returns the radius of the planet in kilometers.
/// The mass passed in is in units of solar masses.
/// This formula is listed as eq.9 in (Fogg 1985), although some typos
/// crop up in that eq. See "The Internal Constitution of Planets", by
/// Dr. D. S. Kothari, Mon. Not. of the Royal Astronomical Society, vol 96
/// pp.833-843, 1936 for the derivation. Specifically, this is Kothari's
/// eq.23, which appears on page 840.
pub fn kothari_radius(mass: f64, mass_type: MassType, zone: OrbitalZone) -> f64 {
    let (atomic_weight, atomic_num) = match zone {
        OrbitalZone::Zone1 => {
            if mass_type == MassType::GasGiant {
                (9.5, 4.5)
            } else {
                (15.0, 8.0)
            }
        }
        OrbitalZone::Zone2 => {
            if mass_type == MassType::GasGiant {
                (2.47, 2.0)
            } else {
                (10.0, 5.0)
            }
        }
        OrbitalZone::Zone3 => {
            if mass_type == MassType::GasGiant {
                (7.0, 4.0)
            } else {
                (10.0, 5.0)
            }
        }
    };

    let mut temp: f64 = atomic_weight * atomic_num;
    temp =
        (2.0 * consts::BETA_20 * consts::SOLAR_MASS_IN_GRAMS.powf(1.0 / 3.0)) / (consts::A1_20 * temp.powf(1.0 / 3.0));

    let mut temp2 = consts::A2_20 * atomic_weight.powf(4.0 / 3.0) * consts::SOLAR_MASS_IN_GRAMS.powf(2.0 / 3.0);
    temp2 *= mass.powf(2.0 / 3.0);
    temp2 /= consts::A1_20 * atomic_num.powf(2.0);
    temp2 += 1.0;
    temp /= temp2;

    (temp * mass.powf(1.0 / 3.0)) / consts::CM_PER_KM
}

/// The mass passed in is in units of solar masses, and the luminosity is
/// a unitless ratio. The density is returned in units of grams/cc.
pub fn empirical_density(mass: f64, orb_radius: f64, mass_type: MassType, star_luminosity: f64) -> f64 {
    let temp = (mass * consts::SUN_MASS_IN_EARTH_MASSES).powf(1.0 / 8.0);
    let temp2 = star_luminosity.sqrt();
    let temp = temp * (temp2 / orb_radius).powf(0.25);
    if mass_type == MassType::GasGiant {
        temp * 1.2
    } else {
        temp * 5.5
    }
}

/// The mass passed in is in units of solar masses, and the equatorial
/// radius is in km. The density is returned in units of grams/cc.
pub fn volume_density(mass: f64, equat_radius: f64) -> f64 {
    let mass_in_grams = mass * consts::SOLAR_MASS_IN_GRAMS;
    let equat_radius_cm = equat_radius * consts::CM_PER_KM;
    let volume = (4.0 * std::f64::consts::PI * equat_radius_cm.powf(3.0)) / 3.0;
    mass_in_grams / volume
}

/// The separation is in units of AU, and both masses are in units of solar
/// masses. The period returned is in terms of Earth days.
pub fn period(separation: f64, small_mass: f64, large_mass: f64) -> f64 {
    let period_in_years = (separation.powf(3.0) / (small_mass + large_mass)).sqrt();
    period_in_years * env::DAYS_IN_A_YEAR
}

/// (Fogg 1985) information for this routine came from (Dole 1964) "Habitable Planets
/// for Man", Blaisdell Publishing Company, NY, 1964. From this, he came
/// up with his eq.12, which is the equation for the `base_angular_velocity`
/// below. He then used an equation for the change in angular velocity per
/// time (dw/dt) from (Goldreich & Soter 1966) "Q in the Solar
/// System" in Icarus, vol 5, pp.375-389 (1966). Using as a comparison the
/// change in angular velocity for the Earth, Fogg has come up with an
/// approximation for our new planet (his eq.13) and take that into account.
/// This is used to find `change_in_angular_velocity` below.
///
/// Input parameters are mass (in solar masses), radius (in Km), orbital
/// period (in days), orbital radius (in AU), density (in g/cc),
/// eccentricity, and whether it is a gas giant or not.
/// The length of the day is returned in units of hours.
///
/// `resonance` is set to `true` if spin resonance occurs.
#[allow(clippy::too_many_arguments)]
pub fn day_length(
    mass: f64,
    radius: f64,
    eccentricity: f64,
    density: f64,
    orb_radius: f64,
    orb_period: f64,
    mass_type: MassType,
    stell_mass_ratio: f64,
    age: f64,
    resonance: &mut bool,
) -> f64 {
    let k2 = if mass_type == MassType::GasGiant { 0.24 } else { 0.33 };
    let planetary_mass_in_grams = mass * consts::SOLAR_MASS_IN_GRAMS;
    let equatorial_radius_in_cm = radius * consts::CM_PER_KM;
    let year_in_hours = orb_period * 24.0;

    let base_angular_velocity =
        (2.0 * env::J * planetary_mass_in_grams / (k2 * equatorial_radius_in_cm.powf(2.0))).sqrt();

    // This next calculation determines how much the planet's rotation is
    // slowed by the presence of the star.
    let change_in_angular_velocity = env::CHANGE_IN_EARTH_ANG_VEL
        * (density / env::EARTH_DENSITY)
        * (equatorial_radius_in_cm / env::EARTH_RADIUS)
        * (env::EARTH_MASS_IN_GRAMS / planetary_mass_in_grams)
        * stell_mass_ratio.powf(2.0)
        * (1.0 / orb_radius.powf(6.0));

    let ang_velocity = base_angular_velocity + (change_in_angular_velocity * age);

    // Now we change from rad/sec to hours/rotation.
    let mut stopped = false;
    let mut day_in_hours = 0.0;
    if ang_velocity <= 0.0 {
        stopped = true;
    } else {
        day_in_hours = env::RADIANS_PER_ROTATION / (env::SECONDS_PER_HOUR * ang_velocity);
    }

    *resonance = false;
    if (day_in_hours >= year_in_hours) || stopped {
        if eccentricity > 0.1 {
            let spin_resonance_factor = (1.0 - eccentricity) / (1.0 + eccentricity);
            *resonance = true;
            return spin_resonance_factor * year_in_hours;
        }
        return year_in_hours;
    }

    day_in_hours
}

/// The orbital radius is expected in units of Astronomical Units (AU).
/// Inclination is returned in units of degrees.
pub fn inclination(orb_radius: f64) -> i32 {
    let temp = (orb_radius.powf(0.2) * about(env::EARTH_AXIAL_TILT, 0.4)) as i32;
    temp % 360
}

/// This function implements the escape velocity calculation. Note that
/// it appears that Fogg's eq.15 is incorrect.
/// The mass is in units of solar mass, the radius in kilometers, and the
/// velocity returned is in cm/sec.
pub fn escape_vel(mass: f64, radius: f64) -> f64 {
    let mass_in_grams = mass * consts::SOLAR_MASS_IN_GRAMS;
    let radius_in_cm = radius * consts::CM_PER_KM;
    (2.0 * env::GRAV_CONSTANT * mass_in_grams / radius_in_cm).sqrt()
}

/// This is Fogg's eq.16. The molecular weight (usually assumed to be N2)
/// is used as the basis of the Root Mean Square (RMS) velocity of the
/// molecule or atom. The velocity returned is in cm/sec.
pub fn rms_vel(molecular_weight: f64, orb_radius: f64, luminosity: f64) -> f64 {
    let exospheric_temp = env::EARTH_EXOSPHERE_TEMP * (luminosity / orb_radius.powf(2.0));
    ((3.0 * env::MOLAR_GAS_CONST * exospheric_temp) / molecular_weight).sqrt() * env::CM_PER_METER
}

/// This function returns the smallest molecular weight retained by the
/// body, which is useful for determining the atmosphere composition.
/// Orbital radius is in A.U.(ie: in units of the earth's orbital radius),
/// mass is in units of solar masses, and equatorial radius is in units of
/// kilometers.
pub fn molecule_limit(mass: f64, equat_radius: f64) -> f64 {
    let esc_velocity = escape_vel(mass, equat_radius);
    (3.0 * (env::GAS_RETENTION_THRESHOLD * env::CM_PER_METER).powf(2.0)
        * env::MOLAR_GAS_CONST
        * env::EARTH_EXOSPHERE_TEMP)
        / esc_velocity.powf(2.0)
}

/// This function calculates the surface acceleration of a planet. The
/// mass is in units of solar masses, the radius in terms of km, and the
/// acceleration is returned in units of cm/sec2.
pub fn accel(mass: f64, radius: f64) -> f64 {
    env::GRAV_CONSTANT * (mass * consts::SOLAR_MASS_IN_GRAMS) / (radius * consts::CM_PER_KM).powf(2.0)
}

/// This function calculates the surface gravity of a planet. The
/// acceleration is in units of cm/sec2, and the gravity is returned in
/// units of Earth gravities.
pub fn gravity(acceleration: f64) -> f64 {
    acceleration / env::EARTH_ACCELERATION
}

/// Note that if the orbital radius of the planet is greater than or equal
/// to R_inner, 99% of it's volatiles are assumed to have been deposited in
/// surface reservoirs (otherwise, it suffers from the greenhouse effect).
pub fn grnhouse(zone: OrbitalZone, orb_radius: f64, r_greenhouse: f64) -> bool {
    matches!(zone, OrbitalZone::Zone1) && (orb_radius < r_greenhouse)
}

/// This implements Fogg's eq.17. The `inventory` returned is unitless.
/// Returns a measure of the amount of gasses locked up inside the planet.
///
/// The zone 3 proportion constant has been raised from 250 to 2500 to
/// better reflect volatile delivery by cometary bombardment in the outer
/// system (the original value produced near-zero hydrospheres on all outer
/// rocky bodies regardless of temperature). The non-greenhouse retention
/// factor is now graduated by zone rather than a flat /100 divisor: inner
/// rocky worlds plausibly retain more volatiles than outer ones without a
/// full runaway greenhouse.
pub fn vol_inventory(
    mass: f64,
    esc_velocity: f64,
    rms_velocity: f64,
    stellar_mass: f64,
    zone: OrbitalZone,
    greenhouse_effect: bool,
) -> f64 {
    vol_inventory_with_profile(
        mass,
        esc_velocity,
        rms_velocity,
        stellar_mass,
        legacy_volatile_retention_profile(zone),
        greenhouse_effect,
    )
}

/// This implements Fogg's eq.18, although it has been changed somewhat
/// to account for planets so close to a star that their atmosphere has
/// been blown off.
///
/// Input                   Units       Description
/// _____                   _____       ___________
/// volatile_gas_inventory  unitless    A measure of the gasses available
/// equat_radius            kilometers  Equatorial radius of the planet
/// grav                    gees        Surface gravity
///
/// Output                  Units       Description
/// _____                   _____       ___________
/// pressure                millibars   Surface atmospheric pressure
pub fn pressure(volatile_gas_inventory: f64, equat_radius: f64, grav: f64) -> f64 {
    let equat_radius = env::KM_EARTH_RADIUS / equat_radius;
    volatile_gas_inventory * grav / equat_radius.powf(2.0)
}

/// This function returns the boiling point of water in an atmosphere of
/// pressure `surf_pressure`, given in millibars. The boiling point is
/// returned in units of Kelvin. This is Fogg's eq.21.
pub fn boiling_point(surf_pressure: f64) -> f64 {
    if surf_pressure <= 0.0 {
        return 0.0;
    }

    let surface_pressure_in_bars = surf_pressure / env::MILLIBARS_PER_BAR;
    1.0 / (surface_pressure_in_bars.ln() / -5050.5 + 1.0 / 373.0)
}

/// This function is Fogg's eq.22. Given the volatile gas inventory and
/// planetary radius of a planet (in Km), this function returns the
/// fraction of the planet covered with water.
/// I have changed the function very slightly: the fraction of Earth's
/// surface covered by water is 71%, not 75% as Fogg used.
pub fn hydro_fraction(volatile_gas_inventory: f64, planet_radius: f64) -> f64 {
    let temp = (0.71 * volatile_gas_inventory / 1000.0) * (env::KM_EARTH_RADIUS / planet_radius).powf(2.0);
    if temp >= 1.0 {
        1.0
    } else {
        temp
    }
}

/// Given the surface temperature of a planet (in Kelvin), this function
/// returns the fraction of cloud cover available. This is Fogg's eq.23.
/// See Hart in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation.
/// This equation is Hart's eq.3.
/// I have modified it slightly using constants and relationships from
/// Glass's book "Introduction to Planetary Geology", p.46.
/// The `CLOUD_COVERAGE_FACTOR` is the amount of surface area on Earth
/// covered by one Kg. of cloud.
pub fn cloud_fraction(surf_temp: f64, smallest_mw_retained: f64, equat_radius: f64, hyd_fraction: f64) -> f64 {
    if smallest_mw_retained > env::WATER_VAPOR {
        return 0.0;
    }

    let surf_area = 4.0 * std::f64::consts::PI * equat_radius.powf(2.0);
    let hydro_mass = hyd_fraction * surf_area * env::EARTH_WATER_MASS_PER_AREA;
    let water_vapor_in_kg = (0.00000001 * hydro_mass) * f64::exp(env::Q2_36 * (surf_temp - 288.0));
    let fraction = env::CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surf_area;

    if fraction >= 1.0 {
        1.0
    } else {
        fraction
    }
}

/// Given the surface temperature of a planet (in Kelvin), this function
/// returns the fraction of the planet's surface covered by ice. This is
/// Fogg's eq.24. See Hart[24] in Icarus vol.33, p.28 for an explanation.
/// I have changed a constant from 70 to 90 in order to bring it more in
/// line with the fraction of the Earth's surface covered with ice, which
/// is approximatly .016 (=1.6%).
pub fn ice_fraction(hyd_fraction: f64, mut surf_temp: f64) -> f64 {
    if surf_temp > 328.0 {
        surf_temp = 328.0;
    }
    let mut temp = ((328.0 - surf_temp) / 90.0).powf(5.0);
    if temp > (1.5 * hyd_fraction) {
        temp = 1.5 * hyd_fraction;
    }
    if temp >= 1.0 {
        1.0
    } else {
        temp
    }
}

/// This is Fogg's eq.19. The Earth-equivalent insolation radius is given in AU,
/// the orbital radius in AU, and the temperature returned is in Kelvin.
pub fn eff_temp(earth_equivalent_insolation_radius_au: f64, orb_radius: f64, albedo: f64) -> f64 {
    (earth_equivalent_insolation_radius_au / orb_radius).sqrt()
        * ((1.0 - albedo) / 0.7).powf(0.25)
        * env::EARTH_EFFECTIVE_TEMP
}

/// This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of
/// Earth's Atmosphere" article. The effective temperature given is in
/// units of Kelvin, as is the rise in temperature produced by the
/// greenhouse effect, which is returned.
pub fn green_rise(optical_depth: f64, effective_temp: f64, surf_pressure: f64) -> f64 {
    let convection_factor =
        env::EARTH_CONVECTION_FACTOR * (surf_pressure / env::EARTH_SURF_PRES_IN_MILLIBARS).powf(0.25);
    ((1.0 + 0.75 * optical_depth).powf(0.25) - 1.0) * effective_temp * convection_factor
}

/// The surface temperature passed in is in units of Kelvin.
/// The cloud adjustment is the fraction of cloud cover obscuring each
/// of the three major components of albedo that lie below the clouds.
pub fn planet_albedo(water_fraction: f64, cld_fraction: f64, ice_frc: f64, surf_pressure: f64) -> f64 {
    let mut rock_fraction = 1.0 - water_fraction - ice_frc;
    let mut components = 0.0;

    if water_fraction > 0.0 {
        components += 1.0;
    }
    if ice_frc > 0.0 {
        components += 1.0;
    }
    if rock_fraction > 0.0 {
        components += 1.0;
    }

    let cloud_adjustment = cld_fraction / components;
    if rock_fraction >= cloud_adjustment {
        rock_fraction -= cloud_adjustment;
    } else {
        rock_fraction = 0.0;
    }

    let mut water_fraction = water_fraction;
    if water_fraction > cloud_adjustment {
        water_fraction -= cloud_adjustment;
    } else {
        water_fraction = 0.0;
    }

    let mut ice_frc = ice_frc;
    if ice_frc > cloud_adjustment {
        ice_frc -= cloud_adjustment;
    } else {
        ice_frc = 0.0;
    }

    let cloud_part = cld_fraction * about(env::CLOUD_ALBEDO, 0.2);
    let rock_part = if surf_pressure == 0.0 {
        rock_fraction * about(env::ROCKY_AIRLESS_ALBEDO, 0.3)
    } else {
        rock_fraction * about(env::ROCKY_ALBEDO, 0.1)
    };
    let water_part = water_fraction * about(env::WATER_ALBEDO, 0.2);
    let ice_part = if surf_pressure == 0.0 {
        ice_frc * about(env::AIRLESS_ICE_ALBEDO, 0.4)
    } else {
        ice_frc * about(env::ICE_ALBEDO, 0.1)
    };

    cloud_part + rock_part + water_part + ice_part
}

/// This function returns the dimensionless quantity of optical depth,
/// which is useful in determining the amount of greenhouse effect on a
/// planet.
pub fn opacity(molecular_weight: f64, surf_pressure: f64) -> f64 {
    let mut optical_depth = 0.0;

    if (0.0..10.0).contains(&molecular_weight) {
        optical_depth += 3.0;
    }
    if (10.0..20.0).contains(&molecular_weight) {
        optical_depth += 2.34;
    }
    if (20.0..30.0).contains(&molecular_weight) {
        optical_depth += 1.0;
    }
    if (30.0..45.0).contains(&molecular_weight) {
        optical_depth += 0.15;
    }
    if (45.0..100.0).contains(&molecular_weight) {
        optical_depth += 0.05;
    }

    if surf_pressure >= (70.0 * env::EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth *= 8.333;
    } else if surf_pressure >= (50.0 * env::EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth *= 6.666;
    } else if surf_pressure >= (30.0 * env::EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth *= 3.333;
    } else if surf_pressure >= (10.0 * env::EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth *= 2.0;
    } else if surf_pressure >= (5.0 * env::EARTH_SURF_PRES_IN_MILLIBARS) {
        optical_depth *= 1.5;
    }

    optical_depth
}

/// The temperature calculated is in degrees Kelvin.
/// Quantities already known which are used in these calculations:
///   planet.molec_weight
///   planet.surf_pressure
///     R_ecosphere
///   planet.a
///   planet.volatile_gas_inventory
///   planet.radius
///   planet.boil_point
/// The `counter` variable used in the iteration loop is used to break
/// out of the loop after 100 iterations - just in case the temperature
/// refuses to converge.
///
/// `planet.hydrosphere` is the total water fraction (liquid + frozen); it is
/// never zeroed by temperature. `planet.ice_cover` holds the frozen portion:
/// the full hydrosphere when the surface is below freezing or below the water
/// triple-point pressure (6.1 mbar), or a partial polar-cap fraction (Fogg
/// eq.24) when above freezing. The liquid fraction used for albedo is derived
/// as `(hydrosphere - ice_cover).max(0.0)`.
pub fn iterate_surface_temp(planet: &mut EnviroProperties, earth_equivalent_insolation_radius_au: f64) {
    // Triple-point pressure of water in millibars. Below this, liquid water
    // cannot exist on the surface regardless of temperature.
    const TRIPLE_POINT_PRESSURE_MB: f64 = 6.1;

    let mut albedo = env::EARTH_ALBEDO;
    let water = hydro_fraction(planet.volatile_gas_inventory, planet.radius_in_km);
    let optical_depth = opacity(planet.molec_weight, planet.surf_pressure);
    let mut eff_water: f64;
    let mut clouds: f64;
    let mut ice: f64;
    let mut new_temp = 0.0;
    let mut counter = 0;

    let (final_effective_temp, final_greenhouse_rise) = loop {
        let mut effective_temp = eff_temp(earth_equivalent_insolation_radius_au, planet.a, albedo);
        let previous_temp = if counter == 0 { effective_temp } else { new_temp };

        let mut greenhs_rise = green_rise(optical_depth, effective_temp, planet.surf_pressure);

        // Watch out for overflow:
        if effective_temp + greenhs_rise > f64::MAX {
            effective_temp = f64::MAX;
            greenhs_rise = 0.0;
            new_temp = f64::MAX;
        } else {
            new_temp = effective_temp + greenhs_rise;
        }

        clouds = cloud_fraction(new_temp, planet.molec_weight, planet.radius_in_km, water);

        // All water is frozen when below the freezing point or when pressure
        // is below the triple point; otherwise use Fogg's eq.24 for partial
        // polar-cap coverage above freezing.
        ice = if new_temp <= env::FREEZING_POINT_OF_WATER || planet.surf_pressure < TRIPLE_POINT_PRESSURE_MB {
            water
        } else {
            ice_fraction(water, new_temp)
        };

        // Liquid fraction is the water budget minus what is frozen.
        eff_water = (water - ice).max(0.0);

        albedo = planet_albedo(eff_water, clouds, ice, planet.surf_pressure);
        counter += 1;

        if (new_temp - previous_temp).abs() <= 1.0 || counter >= env::TEMP_ITERATION_LIMIT {
            break (effective_temp, greenhs_rise);
        }
    };

    planet.hydrosphere = water;
    planet.cloud_cover = clouds;
    planet.ice_cover = ice;
    planet.albedo = albedo;
    planet.effective_temp_k = final_effective_temp.max(0.0);
    planet.greenhouse_rise_k = final_greenhouse_rise.max(0.0);
    planet.surf_temp = new_temp;
}

#[cfg(test)]
mod tests {
    use super::{
        compute_enviro_properties_for_body, greenhouse_effect_for_body, volatile_retention_profile_for_body,
        AtmosphereClass,
    };
    use crate::{
        body::{Body, OrbitalZone},
        condensation::{CondensationFractions, PlanetFormationInputs},
        random::set_rng_seed,
        star::Star,
        types::MassType,
    };
    use std::sync::{Mutex, MutexGuard, OnceLock};

    fn test_rng_guard() -> MutexGuard<'static, ()> {
        static TEST_GUARD: OnceLock<Mutex<()>> = OnceLock::new();
        let mutex = TEST_GUARD.get_or_init(|| Mutex::new(()));
        match mutex.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        }
    }

    fn seed_rng(seed: u64) {
        let applied = set_rng_seed(seed);
        assert_eq!(applied, seed, "test seed should remain deterministic");
    }

    fn formation_inputs(temperature_k: f64, water_ice: f64, volatile_ices: f64, gas: f64) -> PlanetFormationInputs {
        PlanetFormationInputs {
            temperature_k,
            condensation_fractions: CondensationFractions {
                refractory_metal: 0.15,
                silicate_rock: (1.0 - water_ice - volatile_ices - gas - 0.15).max(0.0),
                water_ice,
                volatile_ices,
                gas,
            },
        }
    }

    fn test_star() -> Star {
        Star {
            mass_in_sols: 1.0,
            luminosity_in_sols: 1.0,
            r_ecosphere: 1.0,
            r_greenhouse: 1.5,
            ..Star::default()
        }
    }

    fn test_body(a: f64, mass_in_sols: f64, radius_in_km: f64, formation: Option<PlanetFormationInputs>) -> Body {
        Body {
            a,
            mass_in_sols,
            radius_in_km,
            density_in_grams_per_cc: 5.5,
            mass_type: MassType::Planet,
            orbit_zone: OrbitalZone::Zone2,
            formation_inputs: formation,
            ..Body::default()
        }
    }

    #[test]
    fn volatile_profile_prefers_formation_context_over_compatibility_zone() {
        let dry_world = Body {
            a: 1.0,
            mass_in_sols: 3.0e-6,
            mass_type: MassType::Planet,
            orbit_zone: OrbitalZone::Zone3,
            formation_inputs: Some(formation_inputs(1_100.0, 0.0, 0.0, 0.05)),
            ..Body::default()
        };
        let icy_world = Body {
            a: 1.0,
            mass_in_sols: 3.0e-6,
            mass_type: MassType::Planet,
            orbit_zone: OrbitalZone::Zone1,
            formation_inputs: Some(formation_inputs(40.0, 0.35, 0.25, 0.05)),
            ..Body::default()
        };

        let dry_profile = volatile_retention_profile_for_body(&dry_world);
        let icy_profile = volatile_retention_profile_for_body(&icy_world);

        assert!(dry_profile.greenhouse_candidate);
        assert!(!icy_profile.greenhouse_candidate);
        assert!(icy_profile.proportion_const > dry_profile.proportion_const);
        assert!(icy_profile.retention_divisor < dry_profile.retention_divisor);
    }

    #[test]
    fn greenhouse_candidate_comes_from_climate_class_not_legacy_zone() {
        let star = Star {
            r_greenhouse: 1.5,
            ..Star::default()
        };
        let icy_inner_world = Body {
            a: 1.0,
            mass_in_sols: 3.0e-6,
            mass_type: MassType::Planet,
            orbit_zone: OrbitalZone::Zone1,
            formation_inputs: Some(formation_inputs(40.0, 0.35, 0.25, 0.05)),
            ..Body::default()
        };

        let profile = volatile_retention_profile_for_body(&icy_inner_world);
        assert!(!greenhouse_effect_for_body(&icy_inner_world, &star, profile));
    }

    #[test]
    fn atmosphere_classifies_temperate_rocky_world_as_standard_band() {
        let _guard = test_rng_guard();
        seed_rng(101);
        let star = test_star();
        let body = test_body(1.0, 3.0e-6, 6_371.0, Some(formation_inputs(280.0, 0.08, 0.05, 0.01)));

        let props = compute_enviro_properties_for_body(&body, &star);
        assert!(
            matches!(
                props.atmosphere_class,
                AtmosphereClass::Standard | AtmosphereClass::StandardTainted
            ),
            "temperate rocky world classified as {:?} with {:.2} mb and {:?}",
            props.atmosphere_class,
            props.surf_pressure,
            props.atmosphere_components
        );
        assert!(props.surf_pressure >= 710.0);
        assert!(props.surf_pressure <= 1490.0);
        assert!(props.effective_temp_k.is_finite());
        assert!(props.greenhouse_rise_k.is_finite());
        assert!(props.surf_temp.is_finite());
        assert!((props.effective_temp_k + props.greenhouse_rise_k - props.surf_temp).abs() <= 1.0e-6);
    }

    #[test]
    fn atmosphere_classifies_small_dry_world_as_trace_or_airless() {
        let _guard = test_rng_guard();
        seed_rng(102);
        let star = test_star();
        let body = test_body(0.7, 4.0e-7, 2_600.0, Some(formation_inputs(850.0, 0.0, 0.0, 0.0)));

        let props = compute_enviro_properties_for_body(&body, &star);

        assert!(
            matches!(
                props.atmosphere_class,
                AtmosphereClass::Airless | AtmosphereClass::Trace
            ),
            "small dry world classified as {:?} with {:.2} mb and {:?}",
            props.atmosphere_class,
            props.surf_pressure,
            props.atmosphere_components
        );
        assert!(props.surf_pressure <= 90.0);
    }

    #[test]
    fn atmosphere_classifies_hot_volatile_world_as_steam_or_worse() {
        let _guard = test_rng_guard();
        seed_rng(103);
        let star = Star {
            luminosity_in_sols: 1.4,
            r_ecosphere: 1.18,
            r_greenhouse: 1.8,
            ..test_star()
        };
        let body = test_body(0.45, 6.0e-6, 7_200.0, Some(formation_inputs(420.0, 0.30, 0.25, 0.02)));

        let props = compute_enviro_properties_for_body(&body, &star);
        assert!(
            matches!(
                props.atmosphere_class,
                AtmosphereClass::Steam | AtmosphereClass::Corrosive | AtmosphereClass::Insidious
            ),
            "hot volatile world classified as {:?} with {:.2} mb and {:?}",
            props.atmosphere_class,
            props.surf_pressure,
            props.atmosphere_components
        );
        assert!(props.surf_pressure >= 1500.0);
        assert!(props.cloud_cover >= 0.5);
    }

    #[test]
    fn atmosphere_classifies_gas_envelope_world_as_h2_rich_or_gas_giant() {
        let _guard = test_rng_guard();
        seed_rng(104);
        let star = test_star();
        let body = test_body(1.8, 1.2e-5, 24_000.0, Some(formation_inputs(150.0, 0.04, 0.10, 0.55)));

        let props = compute_enviro_properties_for_body(&body, &star);

        assert!(matches!(
            props.atmosphere_class,
            AtmosphereClass::H2Rich | AtmosphereClass::GasGiant
        ));
        assert!(props.cloud_cover >= 0.5);
    }

    #[test]
    fn atmosphere_classifies_extreme_greenhouse_world_as_corrosive_family() {
        let _guard = test_rng_guard();
        seed_rng(105);
        let star = Star {
            luminosity_in_sols: 2.5,
            r_ecosphere: 1.58,
            r_greenhouse: 2.2,
            ..test_star()
        };
        let body = test_body(0.35, 8.0e-6, 7_500.0, Some(formation_inputs(1_000.0, 0.02, 0.30, 0.01)));

        let props = compute_enviro_properties_for_body(&body, &star);
        assert!(matches!(
            props.atmosphere_class,
            AtmosphereClass::Corrosive | AtmosphereClass::Insidious
        ));
        assert!(props.greenhouse_rise_k > 0.0);
        assert!(props.surf_temp > 500.0);
        assert!(props.surf_pressure >= 1500.0);
    }
}
