use std::fmt;

use crate::{get_log_level, log, random::get_random_number, star::Star};

#[derive(Debug, Clone)]
pub struct StarSystem {
    pub stars: Vec<Star>, // Each system has 1 or more stars
}

// Implement the Display trait for StarSystem
impl fmt::Display for StarSystem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "                         SYSTEM  CHARACTERISTICS\n\n")?;
        for (index, star) in self.stars.iter().enumerate() {
            if index == 0 {
                write!(f, "        PRIMARY STAR\n")?;
                write!(f, "{}", star)?;
                if self.stars.len() > 1 {
                    write!(f, "\nCompanion stars present at:\n")?;
                }
            } else {
                write!(f, "{:2} {} {:7.3} AU\n", index, star.stellar_classification(), star.a)?;
            }
        }

        write!(f, "")
    }
}

// Implement the StarSystem struct
impl StarSystem {
    /// Creates a new `StarSystem` either from a specified star type or randomly generates multiple stars.
    ///
    /// This function generates a `StarSystem` based on the `star_type` input. If `star_type` is empty,
    /// it randomly decides the number of stars (1 to 4) based on predefined probabilities and generates each
    /// randomly with varying orbital radii. If `star_type` is provided, it attempts to create a star system
    /// with a single star matching the specified type. It handles errors in star creation by panicking with
    /// an error message.
    ///
    /// # Parameters
    /// - `star_type`: A string slice that optionally specifies the type of star to create. If empty,
    ///   the function generates a random star system. This string typically matches command-line input.
    /// - `rng_seed`: Reset DETERMINISTIC_RNG with the given seed. If zero, generate a new seed.
    ///
    /// # Panics
    /// - The function panics if it fails to parse the provided `star_type` into a valid star configuration.
    ///
    /// # Examples
    /// ```
    /// // Generate a random star system
    /// let random_system = StarSystem::new("");
    ///
    /// // Generate a specific type of star system
    /// let specific_system = StarSystem::new("G3M/1");
    /// ```
    ///
    /// # Notes
    /// - The random generation of the number of stars mimics real astronomical data suggesting that
    ///   more than half of all stars are members of multiple star systems. The exact distribution used
    ///   here (45% single, 35% binary, 15% trinary, 5% quaternary) is a simplification.
    pub fn new(star_type: &str) -> Self {
        let mut stars: Vec<Star> = Vec::new();
        let log_level = *get_log_level!();

        // If "-t" was used to specify the star type (e.g. G3M/1), generate the requested star. Otherwise randomly generate the system
        if star_type.is_empty() {
            // Create a star system with 1 to 4 stars
            let star_count = match get_random_number(1..=100) {
                1..=45 => 1,
                46..=80 => 2,
                81..=95 => 3,
                _ => 4,
            };

            for index in 0..star_count {
                let orbital_radius_in_au = if index == 0 {
                    0.0
                } else {
                    get_random_number(1.0..=150.0) // 1 - 150 AU
                };

                let star = Star::random(orbital_radius_in_au);
                log!(
                    log_level,
                    1,
                    "Random star: {:.2} solar masses - {}",
                    star.mass_in_sols,
                    star.stellar_classification()
                );
                stars.push(star);
            }
        } else {
            // Create a star system with one star based on the description in the "-t" flag
            let orbital_radius_in_au = 0.0;
            let star;

            match Star::from_str(star_type, orbital_radius_in_au) {
                Ok(value) => star = value,
                Err(err) => panic!("Error: {}", err),
            }

            log!(
                log_level,
                1,
                "Star: {:.2} solar masses - {}",
                star.mass_in_sols,
                star.stellar_classification()
            );
            stars.push(star);
        }

        for star in &mut stars {
            star.accrete();
        }

        StarSystem { stars }
    }
}
