use rand::Rng;
use std::fmt;

use crate::star::Star;

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
                write!(
                    f,
                    "{:2} {} {:7.3} AU\n",
                    index,
                    star.stellar_classification(),
                    star.orbital_radius_in_au
                )?;
            }
        }

        write!(f, "")
    }
}

// Implement the StarSystem struct
impl StarSystem {
    /// Parses command-line input from the "-t" flag
    pub fn from_str(input: &str) -> Result<StarSystem, &'static str> {
        let orbital_radius_in_au = 0.0; // These are never multiple star systems
        match Star::from_str(input, orbital_radius_in_au) {
            Ok(mut star) => {
                star.accrete();
                let stars = vec![star];
                Ok(StarSystem { stars: stars })
            }
            Err(err) => Err(err),
        }
    }

    /// Randomly-generated systems will have 1 - 4 stars.
    /// The percentage of double, triple, and quadruple star systems is basically pulled from a hat - the
    /// best estimates of the actual frequencies of these kind of systems I could find said only that "more
    /// than half of all stars are members of multiple star systems".
    pub fn random() -> Self {
        let star_count;
        match rand::thread_rng().gen_range(1..=100) {
            1..=45 => {
                star_count = 1;
            }
            46..=80 => {
                star_count = 2;
            }
            81..=95 => {
                star_count = 3;
            }
            96..=100 => {
                star_count = 4;
            }
            _ => unreachable!("blah"),
        }

        // Create the stars
        let mut stars: Vec<Star> = Vec::new();
        for index in 0..star_count {
            let orbital_radius_in_au = if index == 0 {
                0.0
            } else {
                rand::thread_rng().gen_range(1.0..=150.0) // 1 - 150 AU
            };

            let mut star = Star::random(orbital_radius_in_au);
            star.accrete();
            stars.push(star);
        }

        StarSystem { stars }
    }
}
