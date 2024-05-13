// src/orbital_body.rs

use crate::types::MassType;
use crate::{consts, get_log_level, log}; // Import the Star struct

#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub a: f64,
    pub e: f64,
    pub mass: f64,
    pub mass_type: MassType,
    pub local_dust_density: f64,
    pub critical_mass_limit: f64, // The mass at which the body begins to accrete gas
}
impl Default for Body {
    fn default() -> Self {
        Body {
            a: 0.0,
            e: 0.0,
            mass: 0.0,
            mass_type: MassType::Planet,
            local_dust_density: 0.0,
            critical_mass_limit: 0.0,
        }
    }
}
impl Body {
    /// Constructs a new `Body` with specified orbital and physical properties.
    ///
    /// This constructor initializes a `Body` with given values for semi-major axis (`a`), eccentricity (`e`),
    /// mass, and type of mass (`mass_type`). The `Body` created represents a celestial body with these characteristics,
    /// which can be a star, planet, or gas giant depending on `mass_type`.
    ///
    /// # Parameters
    /// - `a`: Semi-major axis of the body's orbit in astronomical units (AU). This defines the average distance
    ///        of the body from the central star or point of orbit.
    /// - `e`: Eccentricity of the orbit. A unitless measure that defines the shape of the orbit, where 0 is a
    ///        perfect circle and values approaching 1 indicate more elongated ellipses.
    /// - `mass`: Mass of the body in solar masses.
    /// - `mass_type`: A `MassType` enum indicating whether the body is a star, planet, or gas giant.
    ///
    /// # Returns
    /// Returns a new instance of `Body`.
    ///
    pub fn new(
        a: f64,
        e: f64,
        mass: f64,
        mass_type: MassType,
        local_dust_density: f64,
        critical_mass_limit: f64,
    ) -> Self {
        Body {
            a,
            e,
            mass,
            mass_type,
            local_dust_density,
            critical_mass_limit,
        }
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
    /// - `mass`: The mass in solar masses of the celestial body at the location of interest.
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
    /// ```
    /// let system = StarSystem {
    ///     local_dust_density: 0.1,
    ///     critical_mass_limit: 0.5,
    /// };
    /// let density = system.local_density(0.3);
    /// println!("Local density for mass 0.3 solar masses: {}", density);
    /// ```
    pub fn local_density(&self, mass: f64) -> f64 {
        consts::K * self.local_dust_density / (1.0 + (self.critical_mass_limit / mass).sqrt() * (consts::K - 1.0))
    }

    pub fn is_trivial_mass(&self) -> bool {
        self.mass < consts::TRIVIAL_MASS
    }

    pub fn mass_in_earth_masses(&self) -> f64 {
        self.mass * consts::SUN_MASS_IN_EARTH_MASSES
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
    /// - `a`: Semi-major axis of the orbit in astronomical units (AU).
    /// - `e`: Eccentricity of the orbit, a unitless measure.
    /// - `luminosity`: Luminosity of the central star in solar luminosities.
    ///
    /// # Returns
    /// The critical limit in solar masses, indicating the threshold above which significant gas accretion can occur.
    ///
    pub fn critical_limit(a: f64, e: f64, luminosity: f64) -> f64 {
        let perihelion = a - a * e;
        consts::B * (perihelion * luminosity.sqrt()).powf(-0.75)
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
        let new_a = (self.mass + other.mass) / ((self.mass / self.a) + (other.mass / other.a));

        // Calculate new eccentricity 'e'
        let angular_momentum = self.mass * self.a.sqrt() * (1.0 - self.e.powf(2.0)).sqrt()
            + other.mass * other.a.sqrt() * (1.0 - other.e.powf(2.0)).sqrt();
        let new_angular_momentum = angular_momentum / ((self.mass + other.mass) * new_a.sqrt());
        let new_e_squared = 1.0 - new_angular_momentum.powf(2.0);
        let new_e = if new_e_squared < 0.0 || new_e_squared >= 1.0 {
            0.0
        } else {
            new_e_squared.sqrt()
        };

        // Update the mass
        self.a = new_a;
        self.e = new_e;
        self.mass += other.mass;

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
}
