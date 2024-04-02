use rand::Rng;
use std::fmt;
// pub mod star;
use crate::star::Star;

pub struct StarSystem {
    /// Each system has 1 or more stars
    pub stars: Vec<Star>,
}

impl fmt::Display for StarSystem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"                         SYSTEM  CHARACTERISTICS\n\n")?;
        for (index, star) in self.stars.iter().enumerate() {
            if index == 0 {
                write!(f,"        PRIMARY STAR\n")?;
                write!(f, "{}", star)?;
                if self.stars.len() > 1 {
                    write!(f,"\nCompanion stars present at:\n")?;
                }
            } else {
                write!(f, "{:2} {} {:7.3} AU\n", index, star.stellar_classification(), star.orbital_radius)?;
            }
        }

        write!(f, "")
    }
}

impl StarSystem {
    /// Parses command-line input from the "-t" flag
    pub fn from_str(input: &str) -> Result<StarSystem, &'static str> {
        match Star::from_str(input) {
            Ok(value) => {
                let stars = vec![value];
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

        // Determine basic characteristics of all the stars:
        let mut stars: Vec<Star> = Vec::new();
        for _ in 0..star_count {
            let mut star = Star::random();
            if stars.len() > 0 {
                star.orbital_radius = rand::thread_rng().gen_range(1.0..=150.0);
            }
            stars.push(star);
        };

        StarSystem { stars: stars }
    }
}
