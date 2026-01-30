// src/orbital_body.rs

use crate::accretion_parameters::ACCRETION_PARAMETERS;
use crate::{accretion_disk::AccretionDisk, consts, get_log_level, log, types::MassType};
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
            accretion_disk: None, // Start without an accretion disk
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
        let mut body = Body::default();
        body.a = a;
        body.e = e;
        body.mass_in_sols = mass_in_sols;
        body.mass_type = mass_type;
        body.radius_in_km = radius_in_km;
        body.local_dust_density = body.dust_density(central_mass_in_sols);
        body.critical_mass_limit = body.critical_limit(stellar_luminosity_in_sols);
        body.orbit_zone = body.calculate_orbit_zone(stellar_luminosity_in_sols);
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
        let mass_influence = (mass / (1.0 + mass)).powf(0.25);
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
        consts::B * (perihelion * stellar_luminosity_in_sols.sqrt()).powf(-0.75)
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
        params.dust_density_coefficient * stellar_mass.sqrt() * f64::exp(-consts::ALPHA * self.a.powf(1.0 / consts::N))
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
    /// Given that ρ_M ≈ ρ_m for our purposes, the formula simplifies to approximately:
    /// d = 2.44 * R
    ///
    /// # Parameters:
    /// - `primary_radius_in_au`: The radius of the primary body in the same units as the desired Roche limit.
    ///
    /// # Returns:
    /// - The Roche limit in the same units as the input radius.
    pub fn roche_limit_in_au(&self) -> f64 {
        2.44 * self.radius_in_km / consts::KM_PER_AU
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
        let new_a =
            (self.mass_in_sols + other.mass_in_sols) / ((self.mass_in_sols / self.a) + (other.mass_in_sols / other.a));

        // Calculate new eccentricity 'e'
        let angular_momentum = self.mass_in_sols * self.a.sqrt() * (1.0 - self.e.powf(2.0)).sqrt()
            + other.mass_in_sols * other.a.sqrt() * (1.0 - other.e.powf(2.0)).sqrt();
        let new_angular_momentum = angular_momentum / ((self.mass_in_sols + other.mass_in_sols) * new_a.sqrt());
        let new_e_squared = 1.0 - new_angular_momentum.powf(2.0);
        let new_e = if new_e_squared < 0.0 || new_e_squared >= 1.0 {
            0.0
        } else {
            new_e_squared.sqrt()
        };

        // Update the mass
        self.a = new_a;
        self.e = new_e;
        self.mass_in_sols += other.mass_in_sols;

        // If the protoplanet had the misfortune to collide with a star, update the corresponding star
        if self.mass_type == MassType::Star {
            // TODO: Implement this
            // node.star_ptr.orbit_radius = node.a;
            // node.star_ptr.stell_mass_ratio = node.mass;
        }

        log!(
            *get_log_level!(),
            1,
            "Collision with a {}! ({:.2}, {:.2} -> {:.2})",
            self.mass_type,
            other.a,
            self.a,
            new_a
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

    /// Calculates the radius of a spherical object given its mass and density.
    ///
    /// # Returns
    /// The radius of the object in kilometers.
    ///
    fn calculate_radius_from_density(&self) -> f64 {
        // Convert mass from solar masses to grams
        let mass_in_grams = self.mass_in_sols * consts::SOLAR_MASS_IN_GRAMS;

        // Calculate volume in cubic centimeters (cm³)
        let volume_cm3 = mass_in_grams / self.density_in_grams_per_cc;

        // Calculate radius in centimeters using the formula for the volume of a sphere:
        // volume = (4/3) * π * radius³
        // Solving for radius:
        // radius = ((3 * volume) / (4 * π))^(1/3)
        let radius_cm = ((3.0 * volume_cm3) / (4.0 * std::f64::consts::PI)).powf(1.0 / 3.0);

        // Convert radius from centimeters to kilometers
        let radius_km = radius_cm / consts::CM_PER_KM;

        // Return the radius in kilometers
        radius_km
    }

    /// Calculates the radius of a planet in kilometers using Kothari's formula.
    ///
    /// The mass passed in is in units of solar masses. This formula is based on
    /// Kothari's equation from "The Internal Constitution of Planets" by Dr. D. S. Kothari,
    /// Mon. Not. of the Royal Astronomical Society, vol 96 pp.833-843, 1936.
    /// Specifically, this is Kothari's eq.23, which appears on page 840.
    ///
    /// # Returns
    /// The radius of the planet in kilometers.
    fn calculate_kothari_radius(&self) -> f64 {
        // Determine atomic weight and atomic number based on zone and mass type
        let (atomic_weight, atomic_num): (f64, f64) = match self.orbit_zone {
            OrbitalZone::Zone1 => {
                if self.mass_type == MassType::GasGiant {
                    (9.5, 4.5)
                } else {
                    (15.0, 8.0)
                }
            }
            OrbitalZone::Zone2 => {
                if self.mass_type == MassType::GasGiant {
                    (2.47, 2.0)
                } else {
                    (10.0, 5.0)
                }
            }
            _ => {
                if self.mass_type == MassType::GasGiant {
                    (7.0, 4.0)
                } else {
                    (10.0, 5.0)
                }
            }
        };

        // Calculate temp
        let temp = atomic_weight * atomic_num;
        let temp = (2.0 * consts::BETA_20 * consts::SOLAR_MASS_IN_GRAMS.powf(1.0 / 3.0))
            / (consts::A1_20 * temp.powf(1.0 / 3.0));

        // Calculate temp2
        let mut temp2 = consts::A2_20 * atomic_weight.powf(4.0 / 3.0) * consts::SOLAR_MASS_IN_GRAMS.powf(2.0 / 3.0);
        temp2 *= self.mass_in_sols.powf(2.0 / 3.0);
        temp2 /= consts::A1_20 * atomic_num.powf(2.0);
        temp2 = 1.0 + temp2;

        // Final calculation of temp
        let temp = temp / temp2;
        let temp = (temp * self.mass_in_sols.powf(1.0 / 3.0)) / consts::CM_PER_KM;

        // Return the radius in kilometers
        temp
    }

    /// Calculates the density of a planetary body based on its properties and the luminosity of its star.
    ///
    /// # Arguments
    /// - `luminosity_in_sols`: The luminosity of the star (in solar units).
    ///
    /// # Returns
    /// - `f64`: The calculated density of the planet in units of grams/cc
    fn calculate_empirical_density(&self, luminosity_in_sols: f64) -> f64 {
        let temp = self.mass_in_earth_masses().powf(1.0 / 8.0);
        let temp2 = luminosity_in_sols.sqrt();
        let temp = temp * (temp2 / self.a).powf(0.25);

        if self.mass_type == MassType::GasGiant {
            temp * 1.2
        } else {
            temp * 5.5
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
        let volume_in_cc = (4.0 * std::f64::consts::PI * radius_in_cm.powf(3.0)) / 3.0;

        // Return the density
        mass_in_grams / volume_in_cc
    }

    /// Initializes the planetary object's properties based on its mass type and the luminosity of its star.
    ///
    /// # Arguments
    /// - `luminosity_in_sols`: The luminosity of the star in solar units.
    ///
    /// # Description
    /// This function sets up key properties of a planetary object:
    /// 1. Determines the orbit zone of the planet based on the star's luminosity.
    /// 2. Calculates the planet's density and radius differently depending on whether the planet is a gas giant or not:
    ///    - **Gas Giant**:
    ///      - Uses empirical density calculations to estimate the density in grams per cubic centimeter.
    ///      - Computes the radius in kilometers based on the calculated density.
    ///    - **Non-Gas Giant**:
    ///      - Determines the radius using the Kothari equation.
    ///      - Derives the density based on the calculated volume.
    ///
    /// # Modifies
    /// - `self.orbit_zone`: Sets the orbital zone classification of the planet.
    /// - `self.density_in_grams_per_cc`: Updates the density of the planet.
    /// - `self.radius_in_km`: Updates the radius of the planet.
    ///
    /// # Notes
    /// - Ensure that the planetary object has its `mass_type` field correctly set before calling this function.
    pub fn initialize(&mut self, luminosity_in_sols: f64) {
        self.orbit_zone = self.calculate_orbit_zone(luminosity_in_sols);

        if self.mass_type == MassType::GasGiant {
            self.density_in_grams_per_cc = self.calculate_empirical_density(luminosity_in_sols);
            self.radius_in_km = self.calculate_radius_from_density();
        } else {
            self.radius_in_km = self.calculate_kothari_radius();
            self.density_in_grams_per_cc = self.calculate_density_from_volume();
        }
    }
}
