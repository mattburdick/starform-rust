//! Planetary environment routines for determining things like planet radius, density,
//! mass, surface temp, etc.
//!
//! Originally ported from enviro.c (1991), now partially modernized with
//! class-based atmosphere and climate approximations.

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
}

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

#[derive(Debug, Clone, Copy)]
struct VolatileRetentionProfile {
    proportion_const: f64,
    retention_divisor: f64,
    greenhouse_candidate: bool,
}

#[derive(Debug, Clone, Copy)]
struct AtmosphereSignals {
    eq_temp_k: f64,
    smallest_mw_retained: f64,
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
}

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

fn greenhouse_effect_for_body(body: &Body, star: &Star, profile: VolatileRetentionProfile) -> bool {
    profile.greenhouse_candidate && body.a < star.r_greenhouse
}

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

fn seed_albedo_for_body(body: &Body) -> f64 {
    match body.climate_class() {
        Some(BodyClimateClass::DryRocky) => 0.16,
        Some(BodyClimateClass::TemperateRocky) => 0.28,
        Some(BodyClimateClass::IceRich) => 0.45,
        Some(BodyClimateClass::GasEnvelope) => 0.42,
        None => 0.30,
    }
}

fn build_atmosphere_signals(
    body: &Body,
    star: &Star,
    radius_km: f64,
    smallest_mw_retained: f64,
    volatile_gas_inventory: f64,
    surface_grav_gees: f64,
    greenhouse_effect: bool,
) -> AtmosphereSignals {
    let (water_inventory, volatile_inventory_fraction, gas_fraction) = inferred_composition_fractions(body);
    let climate_class = body.climate_class();
    let eq_temp_k = eff_temp(
        earth_equivalent_insolation_radius_au(star),
        body.a,
        seed_albedo_for_body(body),
    );
    let escape_velocity = escape_vel(body.mass_in_sols, radius_km);
    let rms_velocity = rms_vel(env::MOL_NITROGEN, body.a, star.luminosity_in_sols);
    let retention_strength = if rms_velocity > 0.0 {
        (escape_velocity / rms_velocity / env::GAS_RETENTION_THRESHOLD).clamp(0.2, 5.0)
    } else {
        0.2
    };
    let base_pressure = pressure(volatile_gas_inventory, radius_km, surface_grav_gees);
    let compressed_base_pressure = if base_pressure > 0.0 {
        env::EARTH_SURF_PRES_IN_MILLIBARS * (base_pressure / env::EARTH_SURF_PRES_IN_MILLIBARS).ln_1p() / 4.0
    } else {
        0.0
    };
    let thermal_pressure_floor =
        if eq_temp_k >= 360.0 && (water_inventory >= 0.10 || volatile_inventory_fraction >= 0.15) {
            let reservoir_strength = (0.7 * water_inventory + volatile_inventory_fraction).clamp(0.0, 1.25);
            let heat_term = ((eq_temp_k - 330.0) / 220.0).clamp(0.0, 2.0);
            env::EARTH_SURF_PRES_IN_MILLIBARS
                * reservoir_strength
                * (1.2 + 1.4 * heat_term + if greenhouse_effect { 0.8 } else { 0.0 })
        } else {
            0.0
        };
    let retained_pressure = compressed_base_pressure
        * (0.85 + 0.15 * surface_grav_gees.clamp(0.1, 4.0))
        * (0.75 + 0.25 * retention_strength)
        * (0.90 + 0.80 * volatile_inventory_fraction + 2.00 * gas_fraction);
    let pressure_proxy_mb = retained_pressure.max(thermal_pressure_floor);

    AtmosphereSignals {
        eq_temp_k,
        smallest_mw_retained,
        volatile_gas_inventory,
        surface_grav_gees,
        pressure_proxy_mb,
        mass_earth: body.mass_in_earth_masses(),
        gas_fraction,
        water_inventory,
        volatile_inventory_fraction,
        climate_class,
        greenhouse_effect,
        mass_type: body.mass_type,
    }
}

fn classify_atmosphere(signals: AtmosphereSignals) -> AtmosphereClass {
    let pressure_atm_proxy = signals.pressure_proxy_mb / env::EARTH_SURF_PRES_IN_MILLIBARS;
    let breathable_mw = (18.0..=42.0).contains(&signals.smallest_mw_retained);
    let temperate_world = matches!(signals.climate_class, Some(BodyClimateClass::TemperateRocky));
    let breathable_band = breathable_mw || temperate_world;
    let volatile_rich = signals.volatile_inventory_fraction >= 0.12 || signals.volatile_gas_inventory >= 700.0;
    let hot_volatile =
        (signals.water_inventory >= 0.10 || signals.volatile_inventory_fraction >= 0.15) && signals.eq_temp_k >= 360.0;
    let gas_rich = matches!(signals.climate_class, Some(BodyClimateClass::GasEnvelope)) || signals.gas_fraction >= 0.18;
    let tainted_candidate = matches!(signals.climate_class, Some(BodyClimateClass::DryRocky))
        || signals.eq_temp_k > 325.0
        || signals.eq_temp_k < 235.0
        || (!breathable_mw && !temperate_world)
        || signals.greenhouse_effect;

    if signals.mass_type == MassType::GasGiant || signals.mass_earth >= 30.0 || (gas_rich && pressure_atm_proxy >= 8.0)
    {
        return AtmosphereClass::GasGiant;
    }

    if (signals.surface_grav_gees > 2.0 && pressure_atm_proxy < 0.12)
        || (signals.surface_grav_gees < 0.25 && pressure_atm_proxy > 4.0)
    {
        return AtmosphereClass::Unusual;
    }

    if gas_rich && signals.mass_earth >= 2.0 && pressure_atm_proxy >= 2.5 {
        return AtmosphereClass::H2Rich;
    }

    if signals.eq_temp_k >= 850.0 && (pressure_atm_proxy >= 2.0 || (hot_volatile && signals.greenhouse_effect)) {
        return AtmosphereClass::Insidious;
    }

    if signals.eq_temp_k >= 500.0 && (pressure_atm_proxy >= 1.0 || (hot_volatile && signals.greenhouse_effect)) {
        return AtmosphereClass::Corrosive;
    }

    if hot_volatile
        && (volatile_rich || signals.water_inventory >= 0.20)
        && (pressure_atm_proxy >= 0.5 || signals.greenhouse_effect)
    {
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
        && (signals.smallest_mw_retained < 4.0 || signals.smallest_mw_retained > 80.0 || signals.eq_temp_k < 220.0)
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

fn atmosphere_pressure_multiplier(class: AtmosphereClass) -> f64 {
    match class {
        AtmosphereClass::Airless => 0.0,
        AtmosphereClass::Trace => 0.25,
        AtmosphereClass::VeryThin | AtmosphereClass::VeryThinTainted => 0.60,
        AtmosphereClass::Thin | AtmosphereClass::ThinTainted => 0.85,
        AtmosphereClass::Standard => 1.00,
        AtmosphereClass::StandardTainted => 1.05,
        AtmosphereClass::Dense => 1.25,
        AtmosphereClass::DenseTainted => 1.40,
        AtmosphereClass::Exotic => 1.10,
        AtmosphereClass::Corrosive => 2.50,
        AtmosphereClass::Insidious => 4.50,
        AtmosphereClass::DenseHigh => 1.80,
        AtmosphereClass::ThinLow => 0.45,
        AtmosphereClass::Steam => 2.60,
        AtmosphereClass::H2Rich => 4.50,
        AtmosphereClass::GasGiant => 12.00,
        AtmosphereClass::Unusual => 1.30,
    }
}

fn class_surface_pressure_mb(class: AtmosphereClass, signals: AtmosphereSignals) -> f64 {
    if matches!(class, AtmosphereClass::Airless) {
        return 0.0;
    }

    let (min_mb, max_mb) = atmosphere_pressure_range_mb(class);
    let adjusted_pressure = signals.pressure_proxy_mb
        * atmosphere_pressure_multiplier(class)
        * (0.9 + 0.25 * signals.surface_grav_gees.clamp(0.1, 4.0).sqrt())
        * (1.0
            + signals.gas_fraction
                * if matches!(class, AtmosphereClass::H2Rich | AtmosphereClass::GasGiant) {
                    4.0
                } else {
                    1.5
                })
        * (1.0 + signals.volatile_inventory_fraction * 0.8);

    adjusted_pressure.clamp(min_mb, max_mb)
}

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

/// Computes a planet's environmental properties from its orbital + physical parameters and its primary star.
///
/// This is a convenience wrapper around the legacy `enviro` routines so callers (like the UI) don't need
/// to manually compose the algorithm.
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

    // For the volatile inventory routine, use nitrogen as the baseline molecule.
    let rms_velocity = rms_vel(env::MOL_NITROGEN, body.a, star.luminosity_in_sols);
    let escape_velocity = escape_vel(body.mass_in_sols, radius_km);
    let volatile_gas_inventory = vol_inventory_with_profile(
        body.mass_in_sols,
        escape_velocity,
        rms_velocity,
        star.mass_in_sols,
        climate_profile,
        greenhouse_effect,
    );

    let surface_grav_gees = gravity(accel(body.mass_in_sols, radius_km));
    let signals = build_atmosphere_signals(
        body,
        star,
        radius_km,
        smallest_mw_retained,
        volatile_gas_inventory,
        surface_grav_gees,
        greenhouse_effect,
    );
    let atmosphere_class = classify_atmosphere(signals);
    let surf_pressure = class_surface_pressure_mb(atmosphere_class, signals);
    let boil_point = boiling_point(surf_pressure);

    let mut props = EnviroProperties {
        a: body.a,
        radius_in_km: radius_km,
        molec_weight: smallest_mw_retained,
        surf_pressure,
        volatile_gas_inventory,
        boil_point,
        atmosphere_class,
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
/// This formula is listed as eq.9 in Fogg's article, although some typos
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

/// Fogg's information for this routine came from Dole "Habitable Planets
/// for Man", Blaisdell Publishing Company, NY, 1964. From this, he came
/// up with his eq.12, which is the equation for the `base_angular_velocity`
/// below. He then used an equation for the change in angular velocity per
/// time (dw/dt) from P. Goldreich and S. Soter's paper "Q in the Solar
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
        assert!(matches!(
            props.atmosphere_class,
            AtmosphereClass::Standard | AtmosphereClass::StandardTainted
        ));
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

        assert!(matches!(
            props.atmosphere_class,
            AtmosphereClass::Airless | AtmosphereClass::Trace
        ));
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
        assert!(matches!(
            props.atmosphere_class,
            AtmosphereClass::Steam | AtmosphereClass::Corrosive | AtmosphereClass::Insidious
        ));
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
