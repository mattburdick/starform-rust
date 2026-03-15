// src/orbital_body.rs

use crate::accretion_parameters::ACCRETION_PARAMETERS;
use crate::{
    accretion_disk::AccretionDisk,
    condensation::{CondensationFractions, CondensationThresholds, PlanetFormationInputs},
    consts, get_log_level, log,
    types::MassType,
};
use std::sync::{Arc, RwLock}; // Make sure to import the global reference

#[derive(Debug, Clone)]
pub enum OrbitalZone {
    Zone1,
    Zone2,
    Zone3,
}

#[derive(Debug, Clone)]
pub struct Body {
    pub a: f64,
    pub e: f64,
    pub mass_in_sols: f64,
    pub mass_type: MassType,
    pub radius_in_km: f64,
    pub local_dust_density: f64,
    pub critical_mass_limit: f64, // The mass at which the body begins to accrete gas
    pub orbit_zone: OrbitalZone,
    pub density_in_grams_per_cc: f64,
    pub formation_inputs: Option<PlanetFormationInputs>,
    /// Mass-weighted composition accumulated during dust sweeps and collisions.
    /// Used by coalesce_body to derive the final formation_inputs, capturing
    /// material from the full feeding zone rather than a single orbital radius.
    pub accreted_fractions: Option<CondensationFractions>,
    pub accretion_disk: Option<Arc<RwLock<AccretionDisk>>>,
}
impl Default for Body {
    fn default() -> Self {
        Body {
            a: 0.0,
            e: 0.0,
            mass_in_sols: 0.0,
            mass_type: MassType::Planet,
            radius_in_km: 0.0,
            local_dust_density: 0.0,
            critical_mass_limit: 0.0,
            orbit_zone: OrbitalZone::Zone1,
            density_in_grams_per_cc: 0.0,
            formation_inputs: None,
            accreted_fractions: None,
            accretion_disk: None, // Start without an accretion disk
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BulkPlanetClass {
    IronRichRocky,
    EarthLikeRocky,
    WaterRichWorld,
    SubNeptune,
    GasGiant,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BodyClimateClass {
    DryRocky,
    TemperateRocky,
    IceRich,
    GasEnvelope,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct BulkCompositionSummary {
    metal_fraction: f64,
    rock_fraction: f64,
    ice_fraction: f64,
    gas_fraction: f64,
    bulk_class: BulkPlanetClass,
}

impl BulkCompositionSummary {
    fn from_formation_inputs(
        mass_type: MassType,
        mass_in_earth_masses: f64,
        formation_inputs: PlanetFormationInputs,
    ) -> Self {
        let fractions = formation_inputs.condensation_fractions;
        let solid_total =
            (fractions.refractory_metal + fractions.silicate_rock + fractions.water_ice + fractions.volatile_ices)
                .clamp(0.0, 1.0);
        let gas_fraction = fractions.gas.clamp(0.0, 1.0);

        let (metal_fraction, rock_fraction, ice_fraction) = if solid_total > 0.0 {
            (
                (fractions.refractory_metal / solid_total).clamp(0.0, 1.0),
                (fractions.silicate_rock / solid_total).clamp(0.0, 1.0),
                ((fractions.water_ice + fractions.volatile_ices) / solid_total).clamp(0.0, 1.0),
            )
        } else {
            (0.0, 0.0, 0.0)
        };

        let bulk_class = if mass_type == MassType::GasGiant {
            BulkPlanetClass::GasGiant
        } else if gas_fraction >= 0.5 && mass_in_earth_masses >= 1.5 {
            BulkPlanetClass::SubNeptune
        } else if solid_total <= 0.05 {
            BulkPlanetClass::EarthLikeRocky
        } else if ice_fraction >= 0.45 {
            BulkPlanetClass::WaterRichWorld
        } else if metal_fraction >= 0.30 {
            BulkPlanetClass::IronRichRocky
        } else {
            BulkPlanetClass::EarthLikeRocky
        };

        Self {
            metal_fraction,
            rock_fraction,
            ice_fraction,
            gas_fraction,
            bulk_class,
        }
    }
}

impl Body {
    /// Constructs a new `Body` instance representing a celestial object in the simulation.
    ///
    /// This method initializes a `Body` with specified orbital parameters, physical properties,
    /// and an optional link to an `AccretionDisk`. This function is versatile, allowing for the
    /// creation of various types of celestial bodies such as stars, planets, or gas giants, based
    /// on the passed `mass_type`.
    ///
    /// # Parameters:
    /// - `a`: The semi-major axis of the orbit in astronomical units (AU).
    /// - `e`: The eccentricity of the orbit.
    /// - `mass_in_sols`: The mass of the body in solar masses.
    /// - `mass_type`: The type of the celestial body (e.g., Star, Planet, GasGiant).
    /// - `radius_in_km`: The radius of the body in km.
    /// - `central_mass_in_sols`: The mass of the central object (Star, Planet, GasGiant) in solar masses.
    /// - `stellar_luminosity_in_sols`: The luminosity of the central star in solar luminosities.
    /// - `accretion_disk`: An optional reference-counted, mutable reference to an `AccretionDisk`
    ///   representing the disk in which the body is located.
    ///
    /// # Returns:
    /// A new instance of `Body` fully initialized with the provided values.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a: f64,
        e: f64,
        mass_in_sols: f64,
        mass_type: MassType,
        radius_in_km: f64,
        central_mass_in_sols: f64,
        stellar_luminosity_in_sols: f64,
        accretion_disk: Option<Arc<RwLock<AccretionDisk>>>,
    ) -> Self {
        let mut body = Body {
            a,
            e,
            mass_in_sols,
            mass_type,
            radius_in_km,
            ..Body::default()
        };
        body.local_dust_density = body.dust_density(central_mass_in_sols);
        body.critical_mass_limit = body.critical_limit(stellar_luminosity_in_sols);
        body.orbit_zone = body.calculate_orbit_zone(stellar_luminosity_in_sols);
        if body.mass_type != MassType::Star {
            body.formation_inputs = Some(body.formation_inputs_from_luminosity(stellar_luminosity_in_sols));
        }
        body.accretion_disk = accretion_disk;

        body
    }

    pub fn collects_gas(&self, mass: f64) -> bool {
        mass >= self.critical_mass_limit
    }

    /// Calculates the local density influenced by both dust and gas components based on the provided mass.
    ///
    /// This function determines the density at a location by modifying the base dust density to account for the
    /// presence of gas, especially as the mass of the object at that location approaches or exceeds a critical
    /// mass threshold. The calculation uses an interpolation formula that scales the dust density by a factor
    /// dependent on the mass relative to the critical mass limit.
    ///
    /// # Parameters
    /// - `central_mass_in_sols`: The mass in solar masses of the central body.
    ///
    /// # Returns
    /// Returns the modified local density as a floating-point number, which includes contributions from both dust
    /// and gas, depending on the mass.
    ///
    /// # Formula
    /// The local density \( \rho \) is calculated as:
    /// \[
    /// \rho = \frac{K \cdot \rho_{d}}{1 + \sqrt{\frac{m_{c}}{m}} \cdot (K - 1)}
    /// \]
    /// where:
    /// - `K` is a constant that represents the dust-to-gas ratio.
    /// - `\rho_{d}` is the density of dust at the location.
    /// - `m_{c}` is the critical mass at which significant amounts of gas begin to accumulate.
    ///
    /// # Examples
    /// ```rust
    /// use starform_rust::body::Body;
    ///
    /// let body = Body {
    ///     local_dust_density: 0.1,
    ///     critical_mass_limit: 0.5,
    ///     ..Body::default()
    /// };
    ///
    /// let density = body.local_density(0.3);
    /// println!("Local density: {}", density);
    /// ```
    pub fn local_density(&self, central_mass_in_sols: f64) -> f64 {
        // Lock the global ACCRETION_PARAMETERS and read K
        let params = ACCRETION_PARAMETERS.lock().unwrap();

        let mut ratio_of_gas_to_dust = (1.0 / params.percent_dust_in_cloud) * 100.0;
        ratio_of_gas_to_dust = ratio_of_gas_to_dust.min(1000.0);
        ratio_of_gas_to_dust * self.local_dust_density
            / (1.0 + (self.critical_mass_limit / central_mass_in_sols).sqrt() * (ratio_of_gas_to_dust - 1.0))
    }

    pub fn is_trivial_mass(&self) -> bool {
        self.mass_in_sols < consts::TRIVIAL_MASS
    }

    pub fn mass_in_earth_masses(&self) -> f64 {
        self.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES
    }

    /// Inserts a new `Body` into a sorted vector of `Body` objects, maintaining the order by the semi-major axis `a`.
    ///
    /// This function inserts `body` into the `bodies` vector such that the vector remains sorted in ascending order
    /// based on the `a` attribute of each `Body`. If `body` is less than the `a` of an existing body in the vector,
    /// it is inserted before that body. Otherwise, it is added at the end.
    ///
    /// # Parameters
    /// - `bodies`: A mutable reference to the vector of `Body` objects where the new body will be inserted.
    /// - `body`: The `Body` object to insert into the vector.
    ///
    pub fn insert(bodies: &mut Vec<Body>, body: Body) {
        // Find the first position where 'a' is greater
        let index = bodies.iter().position(|x| body.a < x.a).unwrap_or(bodies.len());

        bodies.insert(index, body);
    }

    pub fn gravitational_effect_limits(a: f64, e: f64, mass: f64) -> (f64, f64, f64, f64) {
        // The protoplanet will capture particles based on its mass relative to the primary star (or planet if a moon).
        let ra = a * (1.0 + e); // aphelion distance (farthest point from the star)
        let rp = a * (1.0 - e); // perihelion distance (closest point to the star)

        // Calculate the distance from the mass where matter is affected by its gravity. This depends on its distance from the star.
        let mass_influence = (mass / (1.0 + mass)).sqrt().sqrt();
        let xa = ra * mass_influence;
        let xp = rp * mass_influence;

        let inner_effect_limit = (rp - xp) / (1.0 + consts::CLOUD_ECCENTRICITY);
        let outer_effect_limit = (ra + xa) / (1.0 - consts::CLOUD_ECCENTRICITY);

        (inner_effect_limit, outer_effect_limit, xp, xa)
    }

    /// Calculates the critical limit of a planetary body based on its orbital parameters and the luminosity of its central star.
    ///
    /// The critical limit refers to the mass at which a planet begins to accrete gas from the protoplanetary disk significantly,
    /// as opposed to just dust. This function uses the perihelion distance (closest approach to the star) to adjust for the effects
    /// of the star's luminosity on the accretion process.
    ///
    /// # Parameters
    /// - `stellar_luminosity_in_sols`: Luminosity of the central star in solar luminosities.
    ///
    /// # Returns
    /// The critical limit in solar masses, indicating the threshold above which significant gas accretion can occur.
    ///
    pub fn critical_limit(&self, stellar_luminosity_in_sols: f64) -> f64 {
        let perihelion = self.a - self.a * self.e;
        let base = perihelion * stellar_luminosity_in_sols.sqrt();
        let base_sqrt = base.sqrt();
        consts::B / (base_sqrt * base_sqrt.sqrt())
    }

    /// Calculates the dust density at a given orbital radius around a star.
    ///
    /// This function computes the dust density based on the stellar mass and the distance from the star,
    /// using an exponential decay model. The formula incorporates predefined constants to adjust the
    /// model based on empirical or theoretical data.
    ///
    /// # Parameters
    /// - `stellar_mass`: The mass of the star in solar masses.
    ///
    /// # Returns
    /// Returns the dust density at the specified orbital radius as a floating-point number.
    ///
    /// # Formula
    /// The dust density \( d \) at a distance \( a \) from a star of mass \( m \) is calculated as:
    /// \[
    /// \rho_{d} = A \cdot \sqrt{M_{\odot}} \cdot e^{-\alpha \cdot a^{1/n}}
    /// \]
    /// where:
    /// - \( A \) is the dust density coefficient (`DUST_DENSITY_COEFF`),
    /// - \( \alpha \) is a decay constant that influences how quickly the dust density decreases with distance (`ALPHA`),
    /// - \( n \) affects the curvature of the decay curve (`N`).
    pub fn dust_density(&self, stellar_mass: f64) -> f64 {
        // Lock the global ACCRETION_PARAMETERS mutex and retrieve the parameters
        let params = ACCRETION_PARAMETERS.lock().unwrap();

        // Use the globally stored parameters
        params.dust_density_coefficient * stellar_mass.sqrt() * f64::exp(-consts::ALPHA * self.a.cbrt())
    }

    /// Calculates the Roche limit of a primary body for a fluid satellite.
    /// The Roche limit is the minimum distance at which a satellite can orbit a primary body without being torn apart by tidal forces.
    /// This function assumes that both the primary and the satellite are of similar density, a condition which simplifies the formula.
    ///
    /// The formula used here is:
    /// d = R * (2 * (ρ_M / ρ_m))^(1/3)
    /// Where:
    /// - R is the radius of the primary body.
    /// - ρ_M and ρ_m are the densities of the primary body and the satellite, respectively.
    ///
    /// Given that ρ_M ≈ ρ_m for our purposes, the formula simplifies to approximately:
    ///
    /// ```text
    /// d = 2.44 * R
    /// ```
    ///
    /// # Parameters:
    /// - `primary_radius_in_au`: The radius of the primary body in the same units as the desired Roche limit.
    ///
    /// # Returns:
    /// - The Roche limit in the same units as the input radius.
    pub fn roche_limit_in_au(&self) -> f64 {
        2.44 * self.radius_in_km / consts::KM_PER_AU
    }

    /// Calculates the Hill radius of this body in astronomical units (AU).
    ///
    /// The Hill sphere is a stability boundary for satellites orbiting a smaller body (e.g. moons orbiting
    /// a planet) in the presence of a much larger primary (e.g. a star). In practice it answers:
    /// "How far from this body can an object orbit and remain gravitationally bound?"
    ///
    /// # Hill sphere formula
    ///
    /// $$r_{Hill} = a \left(\frac{m}{3M}\right)^{1/3}$$
    ///
    /// Where:
    /// - $a$ = this body's semi-major axis around the primary (AU)
    /// - $m$ = this body's mass (solar masses)
    /// - $M$ = the primary body's mass (solar masses)
    ///
    /// # Sphere of Influence (SOI)
    ///
    /// The **Sphere of Influence (SOI)** is calculated using this formula:
    ///
    /// $$r_{SOI} = a \left(\frac{m}{M}\right)^{2/5}$$
    ///
    /// Where:
    /// - **a** = the semi-major axis of the smaller body's orbit around the larger one (i.e. the distance between them)
    /// - **m** = mass of the smaller body (e.g. the planet)
    /// - **M** = mass of the larger body (e.g. the star)
    ///
    /// As a concrete example — Earth's SOI:
    ///
    /// - a = 1 AU (150 million km)
    /// - m = Earth's mass (5.97 × 10²⁴ kg)
    /// - M = Sun's mass (1.99 × 10³⁰ kg)
    ///
    /// The ratio m/M ≈ 3 × 10⁻⁶, and raising that to the 2/5 power gives about 0.006, so:
    ///
    /// **Earth's SOI ≈ 900,000 km** (about 145 Earth radii)
    ///
    /// For comparison, the Moon orbits at ~384,000 km — comfortably inside it.
    ///
    /// # How it differs from the Hill Sphere
    ///
    /// The Hill sphere uses the exponent **1/3** instead of **2/5**. The SOI is often considered more
    /// accurate for practical spacecraft navigation because it models where the smaller body's gravity
    /// dominates the trajectory equations over the primary's perturbations. The Hill sphere is more of
    /// a stability boundary for long-lived satellite orbits.
    pub fn hill_radius_in_au(&self, primary_mass_in_sols: f64) -> f64 {
        fn is_positive(value: f64) -> bool {
            value.partial_cmp(&0.0) == Some(std::cmp::Ordering::Greater)
        }

        if !is_positive(primary_mass_in_sols) || !is_positive(self.mass_in_sols) || !is_positive(self.a) {
            return 0.0;
        }

        self.a * (self.mass_in_sols / (3.0 * primary_mass_in_sols)).cbrt()
    }

    /// Performs a collision between two celestial bodies, updating their orbital and physical properties.
    ///
    /// This function calculates the new orbital semi-major axis (`a`), eccentricity (`e`), and combined mass after
    /// a collision between `self` and another `Body`. The calculations assume perfectly inelastic collisions,
    /// where the total momentum and mass are conserved, but kinetic energy is not necessarily conserved.
    ///
    /// # Parameters
    /// - `other`: A reference to the other `Body` involved in the collision.
    ///
    pub fn collide(&mut self, other: &Body) {
        // Capture pre-collision state for logging
        let self_mass_earth = self.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES;
        let other_mass_earth = other.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES;
        let self_type = self.mass_type;
        let other_type = other.mass_type;
        let self_a = self.a;
        let other_a = other.a;

        let new_a =
            (self.mass_in_sols + other.mass_in_sols) / ((self.mass_in_sols / self.a) + (other.mass_in_sols / other.a));

        // Calculate new eccentricity 'e'
        let angular_momentum = self.mass_in_sols * self.a.sqrt() * (1.0 - self.e.powi(2)).sqrt()
            + other.mass_in_sols * other.a.sqrt() * (1.0 - other.e.powi(2)).sqrt();
        let new_angular_momentum = angular_momentum / ((self.mass_in_sols + other.mass_in_sols) * new_a.sqrt());
        let new_e_squared = 1.0 - new_angular_momentum.powi(2);
        let new_e = if !(0.0..1.0).contains(&new_e_squared) {
            0.0
        } else {
            new_e_squared.sqrt()
        };

        // Merge accreted compositions before the mass update so we can use
        // pre-collision masses as weights (Raymond 2008: radial mixing via
        // giant impacts delivers volatiles to inner planets).
        self.accreted_fractions = match (self.accreted_fractions, other.accreted_fractions) {
            (Some(self_frac), Some(other_frac)) => {
                Some(self_frac.mass_weighted_blend(self.mass_in_sols, &other_frac, other.mass_in_sols))
            }
            (existing @ Some(_), None) | (None, existing @ Some(_)) => existing,
            (None, None) => None,
        };

        // Update the mass
        self.a = new_a;
        self.e = new_e;
        self.mass_in_sols += other.mass_in_sols;
        let resulting_mass_earth = self.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES;

        // If the protoplanet had the misfortune to collide with a star, update the corresponding star
        if self.mass_type == MassType::Star {
            // TODO: Implement this
            // node.star_ptr.orbit_radius = node.a;
            // node.star_ptr.stell_mass_ratio = node.mass;
        }

        log!(
            *get_log_level!(),
            1,
            "Collision: {} ({:.3} Earth masses at {:.3} AU) + {} ({:.3} Earth masses at {:.3} AU) -> surviving {} at {:.3} AU with {:.3} Earth masses",
            self_type,
            self_mass_earth,
            self_a,
            other_type,
            other_mass_earth,
            other_a,
            self.mass_type,
            new_a,
            resulting_mass_earth
        );
    }

    /// Calculates the orbital zone of a planet based on its semi-major axis and the stellar luminosity.
    ///
    /// The orbital zone is determined by comparing the planet's orbital radius to specific boundaries that depend on
    /// the square root of the star's luminosity. The zones are defined as:
    ///
    /// - **Zone 1**: Orbital radius is less than `4 * sqrt(luminosity)`.
    /// - **Zone 2**: Orbital radius is between `4 * sqrt(luminosity)` and `15 * sqrt(luminosity)`.
    /// - **Zone 3**: Orbital radius is greater than or equal to `15 * sqrt(luminosity)`.
    ///
    /// # Parameters
    ///
    /// - `stellar_luminosity_in_sols`: The luminosity of the star relative to the Sun (dimensionless, where the Sun's luminosity is 1.0).
    ///
    /// # Returns
    ///
    /// An `OrbitalZone` enum variant indicating the orbital zone of the planet:
    /// - `OrbitalZone::Zone1`
    /// - `OrbitalZone::Zone2`
    /// - `OrbitalZone::Zone3`
    fn calculate_orbit_zone(&self, stellar_luminosity_in_sols: f64) -> OrbitalZone {
        if self.a < 4.0 * stellar_luminosity_in_sols.sqrt() {
            OrbitalZone::Zone1
        } else if self.a < 15.0 * stellar_luminosity_in_sols.sqrt() {
            OrbitalZone::Zone2
        } else {
            OrbitalZone::Zone3
        }
    }

    fn calculate_compatibility_orbit_zone_from_temperature(temperature_k: f64) -> OrbitalZone {
        let zone1_boundary_temperature_k = 278.0 / 4.0_f64.sqrt();
        let zone2_boundary_temperature_k = 278.0 / 15.0_f64.sqrt();

        if temperature_k > zone1_boundary_temperature_k {
            OrbitalZone::Zone1
        } else if temperature_k > zone2_boundary_temperature_k {
            OrbitalZone::Zone2
        } else {
            OrbitalZone::Zone3
        }
    }

    /// Calculates the density of a spherical object given its mass and radius.
    ///
    /// # Returns
    /// The density of the object in grams/cc.
    fn calculate_density_from_volume(&self) -> f64 {
        // Convert mass to grams
        let mass_in_grams = self.mass_in_sols * consts::SOLAR_MASS_IN_GRAMS;

        // Convert equatorial radius to centimeters
        let radius_in_cm = self.radius_in_km * consts::CM_PER_KM;

        // Calculate volume of the sphere
        let volume_in_cc = (4.0 * std::f64::consts::PI * radius_in_cm.powi(3)) / 3.0;

        // Return the density
        mass_in_grams / volume_in_cc
    }

    fn formation_inputs_from_luminosity(&self, luminosity_in_sols: f64) -> PlanetFormationInputs {
        let orbital_radius_au = self.a.max(0.01);
        let reference_temperature_k = 278.0 * luminosity_in_sols.max(0.01).sqrt().sqrt();
        let temperature_k = reference_temperature_k / orbital_radius_au.sqrt();
        let condensation_fractions =
            CondensationFractions::from_temperature(temperature_k, CondensationThresholds::default(), 25.0);

        PlanetFormationInputs {
            temperature_k,
            condensation_fractions,
        }
    }

    pub(crate) fn climate_class(&self) -> Option<BodyClimateClass> {
        if self.mass_type == MassType::Star {
            return None;
        }

        let formation_inputs = self.formation_inputs?;
        let fractions = formation_inputs.condensation_fractions;
        let ice_fraction = (fractions.water_ice + fractions.volatile_ices).clamp(0.0, 1.0);
        let thresholds = CondensationThresholds::default();
        let gas_rich =
            self.mass_type == MassType::GasGiant || (fractions.gas >= 0.5 && self.mass_in_earth_masses() >= 1.5);

        if gas_rich {
            Some(BodyClimateClass::GasEnvelope)
        } else if ice_fraction >= 0.45 {
            Some(BodyClimateClass::IceRich)
        } else if formation_inputs.temperature_k >= thresholds.water_ice_temp_k && ice_fraction <= 0.10 {
            Some(BodyClimateClass::DryRocky)
        } else {
            Some(BodyClimateClass::TemperateRocky)
        }
    }

    fn radius_earth_radii_to_km(radius_in_earth_radii: f64) -> f64 {
        radius_in_earth_radii * consts::unused_constants::KM_EARTH_RADIUS
    }

    fn calculate_rocky_radius_earth_radii(mass_in_earth_masses: f64, composition_scale: f64) -> f64 {
        let mass = mass_in_earth_masses.max(1.0e-6);
        let base_radius = if mass <= 1.0 {
            mass.powf(0.30)
        } else if mass <= 8.0 {
            mass.powf(0.27)
        } else {
            8.0_f64.powf(0.27) * (mass / 8.0).powf(0.22)
        };

        composition_scale * base_radius
    }

    fn calculate_sub_neptune_radius_earth_radii(mass_in_earth_masses: f64, composition: BulkCompositionSummary) -> f64 {
        let envelope_fraction =
            (composition.gas_fraction * (mass_in_earth_masses / 20.0).clamp(0.15, 1.0) * 0.18).clamp(0.02, 0.18);
        let core_mass = (mass_in_earth_masses * (1.0 - envelope_fraction)).max(0.1);
        let core_radius = Self::calculate_rocky_radius_earth_radii(core_mass, 1.0 + 0.10 * composition.ice_fraction);
        let envelope_boost = 1.15 + 2.2 * envelope_fraction.sqrt();

        (core_radius * envelope_boost).clamp(1.6, 4.8)
    }

    fn calculate_gas_giant_radius_earth_radii(mass_in_earth_masses: f64, composition: BulkCompositionSummary) -> f64 {
        let mass = mass_in_earth_masses.max(10.0);
        let base_radius = if mass <= 20.0 {
            4.0 * (mass / 15.0).powf(0.15)
        } else if mass <= 100.0 {
            4.4 * (mass / 20.0).powf(0.35)
        } else if mass <= 318.0 {
            8.0 * (mass / 100.0).powf(0.28)
        } else {
            11.0 * (mass / 318.0).powf(-0.04)
        };
        let volatile_bias = 0.96 + (0.08 * composition.ice_fraction) + (0.04 * composition.gas_fraction);

        (base_radius * volatile_bias).clamp(3.8, 13.5)
    }

    fn calculate_modern_radius_earth_radii(&self, composition: BulkCompositionSummary) -> f64 {
        let mass_in_earth_masses = self.mass_in_earth_masses().max(1.0e-6);

        match composition.bulk_class {
            BulkPlanetClass::IronRichRocky => {
                Self::calculate_rocky_radius_earth_radii(mass_in_earth_masses, 0.78 + 0.08 * composition.metal_fraction)
            }
            BulkPlanetClass::EarthLikeRocky => {
                Self::calculate_rocky_radius_earth_radii(mass_in_earth_masses, 0.94 + 0.12 * composition.rock_fraction)
            }
            BulkPlanetClass::WaterRichWorld => {
                Self::calculate_rocky_radius_earth_radii(mass_in_earth_masses, 1.05 + 0.15 * composition.ice_fraction)
            }
            BulkPlanetClass::SubNeptune => {
                Self::calculate_sub_neptune_radius_earth_radii(mass_in_earth_masses, composition)
            }
            BulkPlanetClass::GasGiant => {
                Self::calculate_gas_giant_radius_earth_radii(mass_in_earth_masses, composition)
            }
        }
    }

    fn initialize_modern_bulk_properties(&mut self, formation_inputs: PlanetFormationInputs) {
        if self.mass_type == MassType::Star {
            return;
        }

        if self.mass_in_sols <= 0.0 {
            self.radius_in_km = 0.0;
            self.density_in_grams_per_cc = 0.0;
            return;
        }

        let composition = BulkCompositionSummary::from_formation_inputs(
            self.mass_type,
            self.mass_in_earth_masses(),
            formation_inputs,
        );
        let radius_in_earth_radii = self.calculate_modern_radius_earth_radii(composition);

        self.radius_in_km = Self::radius_earth_radii_to_km(radius_in_earth_radii);
        self.density_in_grams_per_cc = self.calculate_density_from_volume();
    }

    pub fn initialize(&mut self, luminosity_in_sols: f64) {
        let formation_inputs = self.formation_inputs_from_luminosity(luminosity_in_sols);
        self.formation_inputs = (self.mass_type != MassType::Star).then_some(formation_inputs);
        self.orbit_zone = self.calculate_orbit_zone(luminosity_in_sols);
        self.initialize_modern_bulk_properties(formation_inputs);
    }

    pub fn initialize_from_formation_inputs(&mut self, formation_inputs: PlanetFormationInputs) {
        self.formation_inputs = (self.mass_type != MassType::Star).then_some(formation_inputs);
        self.orbit_zone = Self::calculate_compatibility_orbit_zone_from_temperature(formation_inputs.temperature_k);
        self.initialize_modern_bulk_properties(formation_inputs);
    }
}

#[cfg(test)]
mod tests {
    use super::{Body, BodyClimateClass, BulkCompositionSummary, BulkPlanetClass, OrbitalZone};
    use crate::{
        condensation::{CondensationFractions, CondensationThresholds, PlanetFormationInputs, RegionThermalProfile},
        consts,
        types::MassType,
    };

    fn assert_close(left: f64, right: f64) {
        assert!((left - right).abs() < 1.0e-12, "left={left} right={right}");
    }

    fn solar_masses_from_earth_masses(mass_in_earth_masses: f64) -> f64 {
        mass_in_earth_masses / consts::SUN_MASS_IN_EARTH_MASSES
    }

    fn formation_inputs_at_temperature(temperature_k: f64) -> PlanetFormationInputs {
        PlanetFormationInputs {
            temperature_k,
            condensation_fractions: CondensationFractions::from_temperature(
                temperature_k,
                CondensationThresholds::default(),
                25.0,
            ),
        }
    }

    #[test]
    fn initialize_and_formation_inputs_share_same_modern_bulk_model() {
        let formation_inputs = RegionThermalProfile::from_host_luminosity(1.0)
            .sample_planet_formation_inputs(1.0, CondensationThresholds::default());
        let template = Body {
            a: 1.0,
            mass_in_sols: 3.0e-6,
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut legacy = template.clone();
        let mut formation_driven = template;

        legacy.initialize(1.0);
        formation_driven.initialize_from_formation_inputs(formation_inputs);

        assert!(matches!(legacy.orbit_zone, OrbitalZone::Zone1));
        assert!(matches!(formation_driven.orbit_zone, OrbitalZone::Zone1));
        assert_close(legacy.radius_in_km, formation_driven.radius_in_km);
        assert_close(legacy.density_in_grams_per_cc, formation_driven.density_in_grams_per_cc);
    }

    #[test]
    fn iron_rich_and_ice_rich_worlds_land_in_different_bulk_classes() {
        let iron_inputs = formation_inputs_at_temperature(1_100.0);
        let ice_inputs = formation_inputs_at_temperature(40.0);
        let template = Body {
            a: 1.0,
            mass_in_sols: solar_masses_from_earth_masses(1.0),
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut iron_world = template.clone();
        let mut ice_world = template;

        let iron_summary = BulkCompositionSummary::from_formation_inputs(
            iron_world.mass_type,
            iron_world.mass_in_earth_masses(),
            iron_inputs,
        );
        let ice_summary = BulkCompositionSummary::from_formation_inputs(
            ice_world.mass_type,
            ice_world.mass_in_earth_masses(),
            ice_inputs,
        );

        iron_world.initialize_from_formation_inputs(iron_inputs);
        ice_world.initialize_from_formation_inputs(ice_inputs);

        assert_eq!(iron_summary.bulk_class, BulkPlanetClass::IronRichRocky);
        assert_eq!(ice_summary.bulk_class, BulkPlanetClass::WaterRichWorld);
        assert!(iron_world.radius_in_km < ice_world.radius_in_km);
        assert!(iron_world.density_in_grams_per_cc > ice_world.density_in_grams_per_cc);
    }

    #[test]
    fn gas_rich_planets_can_expand_into_sub_neptunes() {
        let gas_rich_inputs = formation_inputs_at_temperature(300.0);
        let mostly_solid_inputs = formation_inputs_at_temperature(150.0);
        let template = Body {
            a: 0.8,
            mass_in_sols: solar_masses_from_earth_masses(8.0),
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut sub_neptune = template.clone();
        let mut rocky_analog = template;

        let sub_neptune_summary = BulkCompositionSummary::from_formation_inputs(
            sub_neptune.mass_type,
            sub_neptune.mass_in_earth_masses(),
            gas_rich_inputs,
        );
        let rocky_summary = BulkCompositionSummary::from_formation_inputs(
            rocky_analog.mass_type,
            rocky_analog.mass_in_earth_masses(),
            mostly_solid_inputs,
        );

        sub_neptune.initialize_from_formation_inputs(gas_rich_inputs);
        rocky_analog.initialize_from_formation_inputs(mostly_solid_inputs);

        assert_eq!(sub_neptune_summary.bulk_class, BulkPlanetClass::SubNeptune);
        assert_eq!(rocky_summary.bulk_class, BulkPlanetClass::EarthLikeRocky);
        assert!(sub_neptune.radius_in_km > rocky_analog.radius_in_km);
        assert!(sub_neptune.density_in_grams_per_cc < rocky_analog.density_in_grams_per_cc);
    }

    #[test]
    fn gas_giant_piecewise_radius_remains_category_correct() {
        let cold_outer_inputs = formation_inputs_at_temperature(60.0);
        let neptune_mass = Body {
            a: 20.0,
            mass_in_sols: solar_masses_from_earth_masses(17.0),
            mass_type: MassType::GasGiant,
            ..Body::default()
        };
        let mut neptune_like = neptune_mass.clone();
        let mut jupiter_like = Body {
            a: 5.2,
            mass_in_sols: solar_masses_from_earth_masses(318.0),
            mass_type: MassType::GasGiant,
            ..Body::default()
        };

        neptune_like.initialize_from_formation_inputs(cold_outer_inputs);
        jupiter_like.initialize_from_formation_inputs(cold_outer_inputs);

        assert!(neptune_like.radius_in_km > 20_000.0);
        assert!(neptune_like.radius_in_km < 40_000.0);
        assert!(jupiter_like.radius_in_km > 60_000.0);
        assert!(jupiter_like.radius_in_km < 90_000.0);
        assert!(jupiter_like.radius_in_km > neptune_like.radius_in_km);
        assert!(neptune_like.density_in_grams_per_cc > 0.5);
        assert!(jupiter_like.density_in_grams_per_cc > 0.3);
    }

    #[test]
    fn rocky_radius_grows_monotonically_with_mass() {
        let formation_inputs = formation_inputs_at_temperature(150.0);
        let mut small = Body {
            a: 1.0,
            mass_in_sols: solar_masses_from_earth_masses(0.5),
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut medium = Body {
            a: 1.0,
            mass_in_sols: solar_masses_from_earth_masses(1.0),
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut large = Body {
            a: 1.0,
            mass_in_sols: solar_masses_from_earth_masses(5.0),
            mass_type: MassType::Planet,
            ..Body::default()
        };

        small.initialize_from_formation_inputs(formation_inputs);
        medium.initialize_from_formation_inputs(formation_inputs);
        large.initialize_from_formation_inputs(formation_inputs);

        assert!(small.radius_in_km < medium.radius_in_km);
        assert!(medium.radius_in_km < large.radius_in_km);
        assert!(small.density_in_grams_per_cc > 0.0);
        assert!(large.density_in_grams_per_cc > 0.0);
    }

    #[test]
    fn formation_inputs_temperature_drives_compatibility_zone() {
        let profile = RegionThermalProfile::from_host_luminosity(1.0);
        let mut inner = Body {
            a: 0.5,
            mass_in_sols: 3.0e-6,
            ..Body::default()
        };
        let mut outer = Body {
            a: 20.0,
            mass_in_sols: 3.0e-6,
            ..Body::default()
        };

        inner.initialize_from_formation_inputs(
            profile.sample_planet_formation_inputs(0.5, CondensationThresholds::default()),
        );
        outer.initialize_from_formation_inputs(
            profile.sample_planet_formation_inputs(20.0, CondensationThresholds::default()),
        );

        assert!(matches!(inner.orbit_zone, OrbitalZone::Zone1));
        assert!(matches!(outer.orbit_zone, OrbitalZone::Zone3));
    }

    #[test]
    fn initialization_persists_formation_inputs_and_climate_class() {
        let dry_inputs = formation_inputs_at_temperature(1_100.0);
        let icy_inputs = formation_inputs_at_temperature(40.0);
        let mut dry = Body {
            a: 1.0,
            mass_in_sols: solar_masses_from_earth_masses(1.0),
            mass_type: MassType::Planet,
            ..Body::default()
        };
        let mut icy = dry.clone();

        dry.initialize_from_formation_inputs(dry_inputs);
        icy.initialize_from_formation_inputs(icy_inputs);

        assert_eq!(dry.formation_inputs, Some(dry_inputs));
        assert_eq!(icy.formation_inputs, Some(icy_inputs));
        assert_eq!(dry.climate_class(), Some(BodyClimateClass::DryRocky));
        assert_eq!(icy.climate_class(), Some(BodyClimateClass::IceRich));
    }
}
