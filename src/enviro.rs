//! Planetary environment routines for determining things like planet radius, density,
//! mass, surface temp, etc.
//!
//! Ported from enviro.c (1991) with minimal behavioral changes.

use crate::body::Body;
use crate::body::OrbitalZone;
use crate::consts;
use crate::consts::unused_constants as env;
use crate::random::about;
use crate::star::Star;
use crate::types::MassType;

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
    pub surf_temp: f64,
}

/// Computes a planet's environmental properties from its orbital + physical parameters and its primary star.
///
/// This is a convenience wrapper around the legacy `enviro` routines so callers (like the UI) don't need
/// to manually compose the algorithm.
pub fn compute_enviro_properties_for_body(body: &Body, star: &Star) -> EnviroProperties {
    let zone = body.orbit_zone.clone();

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
    let greenhouse_effect = grnhouse(zone.clone(), body.a, star.r_greenhouse);

    // For the volatile inventory routine, use nitrogen as the baseline molecule.
    let rms_velocity = rms_vel(env::MOL_NITROGEN, body.a, star.luminosity_in_sols);
    let escape_velocity = escape_vel(body.mass_in_sols, radius_km);
    let volatile_gas_inventory = vol_inventory(
        body.mass_in_sols,
        escape_velocity,
        rms_velocity,
        star.mass_in_sols,
        zone.clone(),
        greenhouse_effect,
    );

    let surface_grav_gees = gravity(accel(body.mass_in_sols, radius_km));
    let surf_pressure = pressure(volatile_gas_inventory, radius_km, surface_grav_gees);
    let boil_point = boiling_point(surf_pressure);

    let mut props = EnviroProperties {
        a: body.a,
        radius_in_km: radius_km,
        molec_weight: smallest_mw_retained,
        surf_pressure,
        volatile_gas_inventory,
        boil_point,
        ..Default::default()
    };

    iterate_surface_temp(&mut props, star.r_ecosphere);
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
    let velocity_ratio = esc_velocity / rms_velocity;
    if velocity_ratio < env::GAS_RETENTION_THRESHOLD {
        return 0.0;
    }

    let proportion_const = match zone {
        OrbitalZone::Zone1 => 100000.0,
        OrbitalZone::Zone2 => 75000.0,
        OrbitalZone::Zone3 => 2500.0, // raised from 250; see doc comment above
    };

    let earth_units = mass * consts::SUN_MASS_IN_EARTH_MASSES;
    let temp1 = (proportion_const * earth_units) / stellar_mass;
    let temp2 = about(temp1, 0.2);

    if greenhouse_effect {
        temp2
    } else {
        // Originally a flat /100 for all zones. Now graduated so that inner
        // rocky worlds (Zone1) suffer a smaller penalty than cold outer ones
        // (Zone3), where cold-trapping makes volatile loss more severe.
        let retention_divisor = match zone {
            OrbitalZone::Zone1 => 10.0,
            OrbitalZone::Zone2 => 40.0,
            OrbitalZone::Zone3 => 100.0,
        };
        temp2 / retention_divisor
    }
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

/// This is Fogg's eq.19. The ecosphere radius is given in AU, the orbital
/// radius in AU, and the temperature returned is in Kelvin.
pub fn eff_temp(ecosphere_radius: f64, orb_radius: f64, albedo: f64) -> f64 {
    (ecosphere_radius / orb_radius).sqrt() * ((1.0 - albedo) / 0.7).powf(0.25) * env::EARTH_EFFECTIVE_TEMP
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
pub fn iterate_surface_temp(planet: &mut EnviroProperties, r_ecosphere: f64) {
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

    loop {
        let mut effective_temp = eff_temp(r_ecosphere, planet.a, albedo);
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
            break;
        }
    }

    planet.hydrosphere = water;
    planet.cloud_cover = clouds;
    planet.ice_cover = ice;
    planet.albedo = albedo;
    planet.surf_temp = new_temp;
}
