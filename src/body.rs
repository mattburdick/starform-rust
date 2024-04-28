// src/orbital_body.rs

use crate::{consts, get_log_level}; // Import the Star struct
use crate::types::MassType;

#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub a: f64,
    pub e: f64,
    pub mass: f64,
    pub mass_type: MassType,
}
impl Default for Body {
    fn default() -> Self {
        Body {
            a: 0.0,
            e: 0.0,
            mass: 0.0,
            mass_type: MassType::Planet,
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
    pub fn new(a: f64, e: f64, mass: f64, mass_type: MassType) -> Self {
        Body { a, e, mass, mass_type }
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
        let aphelion = a * (1.0 + e); // rp in Dole's paper
        let perihelion = a * (1.0 - e); // ra in Dole's paper
        let relative_mass_effect = (mass / (1.0 + mass)).powf(0.25);
        let aphelion_influence = aphelion * relative_mass_effect; // xa in Dole's paper
        let parahelion_influence = perihelion * relative_mass_effect; // xp in Dole's paper
        (
            (perihelion - parahelion_influence) / (1.0 + consts::CLOUD_ECCENTRICITY),
            (aphelion + aphelion_influence) / (1.0 - consts::CLOUD_ECCENTRICITY),
            aphelion_influence,
            parahelion_influence,
        )
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
        let new_orbit = (self.mass + other.mass) / ((self.mass / self.a) + (other.mass / other.a));

        // Calculate new eccentricity 'e'
        let angular_momentum = self.mass * self.a.sqrt() * (1.0 - self.e.powf(2.0)).sqrt()
            + other.mass * other.a.sqrt() * (1.0 - other.e.powf(2.0)).sqrt();
        let new_angular_momentum = angular_momentum / ((self.mass + other.mass) * new_orbit.sqrt());
        let new_e_squared = 1.0 - new_angular_momentum.powf(2.0);
        let new_e = if new_e_squared < 0.0 || new_e_squared >= 1.0 {
            0.0
        } else {
            new_e_squared.sqrt()
        };

        // Update the mass
        self.a = new_orbit;
        self.e = new_e;
        self.mass += other.mass;

        // If the protoplanet had the misfortune to collide with a star, update the corresponding star
        if self.mass_type == MassType::Star {
            // TODO: Implement this
            // node.star_ptr.orbit_radius = node.a;
            // node.star_ptr.stell_mass_ratio = node.mass;
        }

        if *get_log_level!() >= 1 {
            let collision_type = match self.mass_type {
                MassType::Star => "star",
                MassType::Planet => "planet",
                MassType::GasGiant => "gas giant",
            };
            println!(
                "  Collision with a {}! ({:.2}, {:.2} -> {:.2})",
                collision_type, other.a, self.a, new_orbit
            );
        }
    }
}
