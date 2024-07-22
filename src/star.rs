// src/star.rs

use lazy_static::lazy_static;
use rand::Rng;
use std::collections::HashMap;
use std::fmt;

use crate::accretion_disk::AccretionDisk;
use crate::body::Body;
use crate::types::MassType;
use crate::{consts, random};

use crate::star::SpectralClass::*;

lazy_static! {
    #[rustfmt::skip]
        static ref STARINFO: HashMap<LuminosityClass, Vec<SpectralInfo>> = {
        let mut star_info = HashMap::new();

        // Main sequence characteristics
        star_info.insert(LuminosityClass::MainSequence, vec![
            SpectralInfo { spec_class: M, spec_num: 9, max_mass: 0.1, percentage: 0.0 },
            SpectralInfo { spec_class: M, spec_num: 5, max_mass: 0.2, percentage: 35.0 },
            SpectralInfo { spec_class: M, spec_num: 0, max_mass: 0.5, percentage: 36.0 },
            SpectralInfo { spec_class: K, spec_num: 5, max_mass: 0.7, percentage: 7.0 },
            SpectralInfo { spec_class: K, spec_num: 0, max_mass: 0.8, percentage: 7.0 },
            SpectralInfo { spec_class: G, spec_num: 5, max_mass: 0.9, percentage: 3.0 },
            SpectralInfo { spec_class: G, spec_num: 0, max_mass: 1.1, percentage: 3.0 },
            SpectralInfo { spec_class: F, spec_num: 5, max_mass: 1.3, percentage: 2.0 },
            SpectralInfo { spec_class: F, spec_num: 0, max_mass: 1.7, percentage: 1.0 },
            SpectralInfo { spec_class: A, spec_num: 5, max_mass: 2.0, percentage: 1.0 },
            SpectralInfo { spec_class: A, spec_num: 0, max_mass: 3.2, percentage: 1.0 },
            SpectralInfo { spec_class: B, spec_num: 5, max_mass: 6.5, percentage: 1.0 },
            SpectralInfo { spec_class: B, spec_num: 0, max_mass: 17.8, percentage: 1.0 },
            SpectralInfo { spec_class: O, spec_num: 5, max_mass: 39.8, percentage: 1.0 },
            SpectralInfo { spec_class: O, spec_num: 0, max_mass: 60.0, percentage: 1.0 },
        ]);

        // White dwarf characteristics
        star_info.insert(LuminosityClass::WhiteDwarf, vec![
            SpectralInfo { spec_class: M, spec_num: 9, max_mass: 0.0, percentage: 0.0 },
            SpectralInfo { spec_class: M, spec_num: 5, max_mass: 0.2, percentage: 0.0 },
            SpectralInfo { spec_class: M, spec_num: 0, max_mass: 0.4, percentage: 0.0 },
            SpectralInfo { spec_class: K, spec_num: 5, max_mass: 0.4, percentage: 1.0 },
            SpectralInfo { spec_class: K, spec_num: 0, max_mass: 0.4, percentage: 1.0 },
            SpectralInfo { spec_class: G, spec_num: 5, max_mass: 0.5, percentage: 1.0 },
            SpectralInfo { spec_class: G, spec_num: 0, max_mass: 0.6, percentage: 1.0 },
            SpectralInfo { spec_class: F, spec_num: 5, max_mass: 0.7, percentage: 4.0 },
            SpectralInfo { spec_class: F, spec_num: 0, max_mass: 0.8, percentage: 8.0 },
            SpectralInfo { spec_class: A, spec_num: 5, max_mass: 1.0, percentage: 28.0 },
            SpectralInfo { spec_class: A, spec_num: 0, max_mass: 0.5, percentage: 32.0 },
            SpectralInfo { spec_class: B, spec_num: 5, max_mass: 0.4, percentage: 13.0 },
            SpectralInfo { spec_class: B, spec_num: 0, max_mass: 0.4, percentage: 9.0 },
            SpectralInfo { spec_class: O, spec_num: 5, max_mass: 0.5, percentage: 1.0 },
            SpectralInfo { spec_class: O, spec_num: 0, max_mass: 0.7, percentage: 1.0 },
            ]);

        // Giant characteristics
        star_info.insert(LuminosityClass::Giant, vec![
            SpectralInfo { spec_class: M, spec_num: 9, max_mass: 8.7, percentage: 0.0 },
            SpectralInfo { spec_class: M, spec_num: 5, max_mass: 7.9, percentage: 12.0 },
            SpectralInfo { spec_class: M, spec_num: 0, max_mass: 6.3, percentage: 19.0 },
            SpectralInfo { spec_class: K, spec_num: 5, max_mass: 5.0, percentage: 26.0 },
            SpectralInfo { spec_class: K, spec_num: 0, max_mass: 4.0, percentage: 25.0 },
            SpectralInfo { spec_class: G, spec_num: 5, max_mass: 3.2, percentage: 5.0 },
            SpectralInfo { spec_class: G, spec_num: 0, max_mass: 2.5, percentage: 4.0 },
            SpectralInfo { spec_class: F, spec_num: 5, max_mass: 2.4, percentage: 2.0 },
            SpectralInfo { spec_class: F, spec_num: 0, max_mass: 2.5, percentage: 1.0 },
            SpectralInfo { spec_class: A, spec_num: 5, max_mass: 2.7, percentage: 1.0 },
            SpectralInfo { spec_class: A, spec_num: 0, max_mass: 3.4, percentage: 1.0 },
            SpectralInfo { spec_class: B, spec_num: 5, max_mass: 7.0, percentage: 1.0 },
            SpectralInfo { spec_class: B, spec_num: 0, max_mass: 30.3, percentage: 1.0 },
            SpectralInfo { spec_class: O, spec_num: 5, max_mass: 60.0, percentage: 1.0 },
            SpectralInfo { spec_class: O, spec_num: 0, max_mass: 70.0, percentage: 1.0 },
            ]);

        // Supergiant characteristics
        star_info.insert(LuminosityClass::Supergiant, vec![
            SpectralInfo { spec_class: M, spec_num: 9, max_mass: 22.3, percentage: 0.0 },
            SpectralInfo { spec_class: M, spec_num: 5, max_mass: 19.9, percentage: 12.0 },
            SpectralInfo { spec_class: M, spec_num: 0, max_mass: 15.8, percentage: 13.0 },
            SpectralInfo { spec_class: K, spec_num: 5, max_mass: 15.0, percentage: 3.0 },
            SpectralInfo { spec_class: K, spec_num: 0, max_mass: 12.6, percentage: 4.0 },
            SpectralInfo { spec_class: G, spec_num: 5, max_mass: 11.6, percentage: 3.0 },
            SpectralInfo { spec_class: G, spec_num: 0, max_mass: 10.0, percentage: 3.0 },
            SpectralInfo { spec_class: F, spec_num: 5, max_mass: 11.8, percentage: 8.0 },
            SpectralInfo { spec_class: F, spec_num: 0, max_mass: 12.6, percentage: 7.0 },
            SpectralInfo { spec_class: A, spec_num: 5, max_mass: 13.2, percentage: 6.0 },
            SpectralInfo { spec_class: A, spec_num: 0, max_mass: 15.8, percentage: 6.0 },
            SpectralInfo { spec_class: B, spec_num: 5, max_mass: 30.2, percentage: 12.0 },
            SpectralInfo { spec_class: B, spec_num: 0, max_mass: 50.1, percentage: 13.0 },
            SpectralInfo { spec_class: O, spec_num: 5, max_mass: 70.0, percentage: 4.0 },
            SpectralInfo { spec_class: O, spec_num: 0, max_mass: 90.0, percentage: 6.0 },
            ]);

        star_info
    };

}
#[derive(Debug, Copy, Clone)]
struct SpectralInfo {
    spec_class: SpectralClass,
    spec_num: i32,
    max_mass: f64,
    percentage: f64,
}

impl SpectralInfo {
    // Get a random entry from the appropriate STARINFO bucket (main sequence, giant, supergiant, etc.)
    // The percentage probability of each entry is used to determine the likelihood of selecting it
    fn get_random(luminosity_class: LuminosityClass) -> SpectralInfo {
        // Generate a random number to select one of the stellar characteristic structs from the appropriate category in STARINFO
        let mut random_percent = rand::thread_rng().gen_range(0.0..=100.0);

        let mut spectral_info: SpectralInfo = SpectralInfo {
            spec_class: G,
            spec_num: 3,
            max_mass: 1.0,
            percentage: 1.0,
        };

        let mut selected_info: Option<SpectralInfo> = None; // Use this to hold the selected item

        // Look up stellar characteristics for the given luminosity class
        if let Some(items) = STARINFO.get(&luminosity_class) {
            for item in items {
                random_percent -= item.percentage;
                if random_percent <= 0.0 {
                    selected_info = Some(*item);
                    break;
                }
            }
        }

        if let Some(info) = selected_info {
            spectral_info = info; // Update the outer variable with the selected item
        }

        spectral_info
    }

    // Get the maximum mass for a given luminosity class, spectral class, and spectral number
    fn get_max_mass(
        luminosity_class: LuminosityClass,
        spec_class: SpectralClass,
        spec_num: i32,
    ) -> Result<f64, &'static str> {
        let mut mass = 0.0;

        if let Some(items) = STARINFO.get(&luminosity_class) {
            for item in items {
                if spec_class == item.spec_class && spec_num >= item.spec_num {
                    mass = item.max_mass;
                    break;
                }
            }
        }

        return if mass == 0.0 {
            Err("No mass found for spectral class")
        } else {
            Ok(mass)
        };
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum SpectralClass {
    O,
    B,
    A,
    F,
    G,
    K,
    M,
    R,
    N,
    S,
}

impl SpectralClass {
    pub fn to_char(&self) -> char {
        match self {
            SpectralClass::O => 'O',
            SpectralClass::B => 'B',
            SpectralClass::A => 'A',
            SpectralClass::F => 'F',
            SpectralClass::G => 'G',
            SpectralClass::K => 'K',
            SpectralClass::M => 'M',
            SpectralClass::R => 'R',
            SpectralClass::N => 'N',
            SpectralClass::S => 'S',
        }
    }
}

impl fmt::Display for SpectralClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SpectralClass::O => write!(f, "O"),
            SpectralClass::B => write!(f, "B"),
            SpectralClass::A => write!(f, "A"),
            SpectralClass::F => write!(f, "F"),
            SpectralClass::G => write!(f, "G"),
            SpectralClass::K => write!(f, "K"),
            SpectralClass::M => write!(f, "M"),
            SpectralClass::R => write!(f, "R"),
            SpectralClass::N => write!(f, "N"),
            SpectralClass::S => write!(f, "S"),
        }
    }
}

impl TryFrom<char> for SpectralClass {
    type Error = &'static str;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'O' => Ok(SpectralClass::O),
            'B' => Ok(SpectralClass::B),
            'A' => Ok(SpectralClass::A),
            'F' => Ok(SpectralClass::F),
            'G' => Ok(SpectralClass::G),
            'K' => Ok(SpectralClass::K),
            'M' => Ok(SpectralClass::M),
            'R' => Ok(SpectralClass::R),
            'N' => Ok(SpectralClass::N),
            'S' => Ok(SpectralClass::S),
            _ => Err("Unknown spectral class"),
        }
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum LuminosityClass {
    Supergiant,
    BrightGiant,
    Giant,
    Subgiant,
    MainSequence,
    WhiteDwarf,
}

impl LuminosityClass {
    // fn to_char(&self) -> char {
    //     match self {
    //         LuminosityClass::Supergiant => 'S',
    //         LuminosityClass::BrightGiant => 'B',
    //         LuminosityClass::Giant => 'G',
    //         LuminosityClass::Subgiant => 'g',
    //         LuminosityClass::MainSequence => 'M',
    //         LuminosityClass::WhiteDwarf => 'D',
    //     }
    // }
    pub fn to_string(&self) -> &'static str {
        match self {
            LuminosityClass::Supergiant => "I (supergiant)",
            LuminosityClass::BrightGiant => "II (bright giant)",
            LuminosityClass::Giant => "III (giant)",
            LuminosityClass::Subgiant => "IV (subgiant)",
            LuminosityClass::MainSequence => "V (main sequence)",
            LuminosityClass::WhiteDwarf => "D (white dwarf)",
        }
    }
}

impl fmt::Display for LuminosityClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.to_string().fmt(f)
    }
}

impl TryFrom<char> for LuminosityClass {
    type Error = &'static str;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'S' => Ok(LuminosityClass::Supergiant),
            'B' => Ok(LuminosityClass::BrightGiant),
            'G' => Ok(LuminosityClass::Giant),
            'g' => Ok(LuminosityClass::Subgiant),
            'M' => Ok(LuminosityClass::MainSequence),
            'D' => Ok(LuminosityClass::WhiteDwarf),
            _ => Err("Unknown luminosity class"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Star {
    pub luminosity_class: LuminosityClass, // Example: MainSequence, Giant, etc.
    pub spectral_class: SpectralClass,     // Example: O, B, A, F, G, K, M
    pub spectral_number: i32,              // Example: 3
    pub temperature_in_kelvin: f64,        // In Kelvin
    pub orbital_radius_in_au: f64,
    pub max_mass_in_sols: f64,
    pub mass_in_sols: f64,
    pub radius_in_au: f64,
    pub luminosity_in_sols: f64,
    pub age: f64,
    pub main_seq_life: f64,
    pub r_ecosphere: f64,
    pub r_greenhouse: f64,
    pub accretion_disk: AccretionDisk,
}

impl fmt::Display for Star {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Stellar Classification:        {}\n", self.stellar_classification())?;
        write!(f, "Stellar mass:                {:>7.3} sols\n", self.mass_in_sols)?;
        write!(
            f,
            "Stellar radius:              {:>7.3} sols ({:>0.4} AU)\n",
            self.radius_in_au,
            self.radius_in_au * consts::SOLAR_RADII_PER_AU
        )?;
        write!(
            f,
            "Stellar luminosity:          {:>7.3} sols \n",
            self.luminosity_in_sols
        )?;
        write!(
            f,
            "Age:                         {:>7.3} billion years\n",
            self.age / 1.0e9
        )?;

        if self.luminosity_class == LuminosityClass::MainSequence {
            write!(
                f,
                "Years left on Main Sequence: {:>7.3} billion years\n",
                (self.main_seq_life - self.age) / 1.0e9
            )?;
        }

        write!(f, "Earthlike insolation at:     {:>7.3} AU\n", self.r_ecosphere)
    }
}

impl Star {
    pub fn stellar_classification(&self) -> String {
        format!(
            "{}{} {}",
            self.spectral_class,
            self.spectral_number,
            self.luminosity_class.to_string()
        )
    }
    /// Sets the spectral and luminosity classifications for a star and updates its properties.
    ///
    /// This method updates the star's spectral class, spectral number, and luminosity class with the given values.
    /// It then recalculates the star's physical properties based on these new classifications.
    ///
    /// # Arguments
    /// * `spectral_class` - The spectral class of the star (e.g., O, B, A, F, G, K, M, R, N, S).
    /// * `spectral_number` - An integer representing the subclass within the spectral class,
    ///   typically ranging from 0 to 9, where 0 is hottest and 9 is coolest within the class.
    /// * `luminosity_class` - The luminosity class of the star (e.g., Main Sequence, Giant,
    ///   Subgiant, Supergiant, Bright Giant, White Dwarf).
    ///
    /// # Example
    /// ```
    /// let mut star = Star::new();
    /// star.set_classifications(SpectralClass::G, 2, LuminosityClass::MainSequence);
    /// ```
    /// After calling this method, the star's `calculate_properties` method is called to update its temperature,
    /// radius, and other related physical properties according to the new classifications.
    pub fn set_classifications(
        &mut self,
        spectral_class: SpectralClass,
        spectral_number: i32,
        luminosity_class: LuminosityClass,
    ) {
        self.spectral_class = spectral_class;
        self.spectral_number = spectral_number;
        self.luminosity_class = luminosity_class;

        self.calculate_properties();
    }

    /// Calculates the luminosity of a star relative to the sun based on its mass and luminosity class.
    ///
    /// This static method implements the mass-luminosity relationship, which is expressed in the
    /// form of a power law as described in equation 3.52 of "Astrophysics I" by Bowers and Deeming.
    /// The `mass_in_sols` is a unitless ratio representing the star's mass relative to the solar mass.
    /// The function adjusts the calculation based on different luminosity classes, using class-specific
    /// constants alpha and beta which are dimensionless.
    ///
    /// # Parameters
    /// * `mass_in_sols` - The star's mass expressed as a multiple of the solar mass.
    /// * `luminosity_class` - The luminosity classification of the star which impacts the mass-luminosity
    ///   relationship. It can be one of the following: MainSequence, Giant, Subgiant, Supergiant,
    ///   BrightGiant, or WhiteDwarf.
    ///
    /// # Returns
    /// * `luminosity_in_sols` - The star's luminosity expressed as a multiple of the solar luminosity.
    ///
    /// # Examples
    /// ```
    /// let luminosity = Star::luminosity_in_sols(1.0, LuminosityClass::MainSequence);
    /// println!("Luminosity of a main-sequence star with solar mass: {}", luminosity);
    /// ```
    ///
    /// # Notes
    /// - For a main-sequence G3 star like the Sun, this function slightly overestimates the luminosity.
    /// - The function generally fits the mass-luminosity curve well across different classes, providing
    ///   an approximation that is useful for various astrophysical calculations and simulations.
    pub fn luminosity_in_sols(mass_in_sols: f64, luminosity_class: LuminosityClass) -> f64 {
        let luminosity_in_sols: f64;

        let log_mass_ratio = mass_in_sols.log10();
        match luminosity_class {
            LuminosityClass::MainSequence => {
                let alpha: f64;
                let beta: f64;

                if mass_in_sols <= 0.5 {
                    alpha = 2.85;
                    beta = -0.15;
                } else if mass_in_sols < 2.5 {
                    alpha = 3.6;
                    beta = 0.073;
                } else {
                    alpha = 2.91;
                    beta = 0.479;
                }
                luminosity_in_sols = 10f64.powf(beta + alpha * log_mass_ratio);
            }
            LuminosityClass::Giant | LuminosityClass::Subgiant => {
                luminosity_in_sols = 10f64.powf(log_mass_ratio * 3.3);
            }
            LuminosityClass::Supergiant | LuminosityClass::BrightGiant => {
                luminosity_in_sols = 10f64.powf((log_mass_ratio + 0.22) / 0.33);
            }
            LuminosityClass::WhiteDwarf => {
                luminosity_in_sols = mass_in_sols * 5.67E-4;
            }
        }

        luminosity_in_sols
    }

    /// Calculates the radius of a star in astronomical units (AU) based on its mass relative to the sun,
    /// luminosity class, and spectral class, following the method described in equation 3.53 of "Astrophysics I"
    /// by Bowers and Deeming.
    ///
    /// This static method uses a logarithmic relationship to compute the star's radius. The formula varies
    /// depending on the luminosity and spectral class, accounting for differences in stellar structure and
    /// evolution.
    ///
    /// # Parameters
    /// * `mass_in_sols` - The mass of the star expressed as a ratio relative to the solar mass.
    /// * `luminosity_class` - The luminosity classification of the star, which significantly affects the radius calculation.
    ///   It can be MainSequence, Giant, Subgiant, Supergiant, BrightGiant, or WhiteDwarf.
    /// * `spectral_class` - The spectral classification of the star, affecting calculations especially for cooler,
    ///   later-type stars such as K and M classes.
    ///
    /// # Returns
    /// * `radius_in_au` - The computed radius of the star in astronomical units, using the solar radius as a conversion factor.
    ///
    /// # Example
    /// ```
    /// let radius = Star::radius_in_au(1.0, LuminosityClass::MainSequence, SpectralClass::G);
    /// println!("Radius of a main-sequence G-type star with solar mass: {}", radius);
    /// ```
    ///
    /// # Notes
    /// - The function applies different scaling factors based on luminosity class and adjusts for cool star types (K, M)
    ///   in the cases of supergiants and bright giants.
    /// - The function for White Dwarfs returns a randomized value around a typical white dwarf radius due to their
    ///   highly compressed nature.
    pub fn radius_in_au(mass_in_sols: f64, luminosity_class: LuminosityClass, spectral_class: SpectralClass) -> f64 {
        let cool_star = spectral_class == K || spectral_class == M;

        let log_mass_ratio = mass_in_sols.log10();
        let radius_in_sols = match luminosity_class {
            LuminosityClass::MainSequence => {
                if mass_in_sols <= 0.4 {
                    10f64.powf(log_mass_ratio + 0.1)
                } else {
                    10f64.powf(0.73 * log_mass_ratio)
                }
            }
            LuminosityClass::Giant | LuminosityClass::Subgiant => 10f64.powf(log_mass_ratio * 2.0),
            LuminosityClass::Supergiant | LuminosityClass::BrightGiant => {
                if cool_star {
                    10f64.powf((log_mass_ratio - 0.32) / 0.34)
                } else {
                    10f64.powf((log_mass_ratio - 2.7) / -0.86)
                }
            }
            LuminosityClass::WhiteDwarf => random::about(0.02, 0.005),
        };

        radius_in_sols * consts::SOLAR_RADII_PER_AU
    }

    /// Calculates the approximate temperature of a star given its luminosity and radius.
    ///
    /// The temperature is estimated based on a modified form of the Stefan-Boltzmann law,
    /// which relates the luminosity of a star to its effective temperature and radius.
    /// This function simplifies the exact relation by using a power law approximation
    /// and scales the result to the temperature of the Sun.
    ///
    /// # Returns
    /// The estimated effective temperature of the star in Kelvin.
    ///
    /// This function uses constants from the module `consts` where `SOLAR_RADII_PER_AU` is defined as
    /// the number of solar radii per astronomical unit and `SOLAR_TEMPERATURE_IN_KELVIN` is the surface
    /// temperature of the sun in Kelvin.
    pub fn temperature_in_kelvin(luminosity_in_sols: f64, radius_in_au: f64) -> f64 {
        let temperature_in_sols = luminosity_in_sols.powf(0.25) / (radius_in_au / consts::SOLAR_RADII_PER_AU).powf(0.5);
        temperature_in_sols * consts::SOLAR_TEMPERATURE_IN_KELVIN
    }

    /// Estimates the main sequence lifetime of a star based on its mass relative to the sun's mass and its
    /// luminosity relative to the sun's luminosity, using a simplified model from stellar astrophysics.
    ///
    /// This function computes the star's main sequence lifetime, assuming that the larger the mass and the higher
    /// the luminosity, the shorter the main sequence phase will be. The formula used is a straightforward
    /// proportionality scaled by a constant representing the sun's estimated main sequence lifetime.
    ///
    /// # Parameters
    /// * `mass_in_sols` - The mass of the star expressed as a ratio relative to the solar mass.
    /// * `luminosity_in_sols` - The luminosity of the star expressed as a ratio relative to the solar luminosity.
    ///
    /// # Returns
    /// * `lifetime` - The estimated main sequence lifetime of the star in years. This is bounded by a minimum value
    ///   to ensure the lifetime never falls below a realistic threshold, reflecting observational constraints.
    ///
    /// # Example
    /// ```
    /// let life_span = Star::main_seq_life(1.0, 1.0);
    /// println!("Main sequence lifetime of a star with solar mass and luminosity: {}", life_span);
    /// ```
    ///
    /// # Notes
    /// - The main sequence lifetime is calculated as a product of a constant (approximately 11 billion years),
    ///   and the ratio of the star's mass to its luminosity. This formula assumes that higher mass and luminosity
    ///   result in faster fuel consumption and thus a shorter lifespan.
    /// - The function ensures that the lifetime does not drop below 10 billion years, which serves as a
    ///   practical lower bound for the main sequence lifetime based on current astrophysical understanding.
    pub fn main_seq_life(mass_in_sols: f64, luminosity_in_sols: f64) -> f64 {
        let lifetime = 1.1E10_f64 * (mass_in_sols / luminosity_in_sols);
        lifetime.max(1.0E10_f64)
    }

    /// Estimates the current age of a star by generating a random value within a range determined by its
    /// calculated main sequence lifetime. This method considers various scenarios based on the longevity of the star.
    ///
    /// The age is generated randomly to simulate the unpredictability of exactly when within its possible
    /// lifespan a star might currently be. This approach is typical in simulations where exact ages are not
    /// determinable but where an age within a realistic range is sufficient.
    ///
    /// # Parameters
    /// * `mass_in_sols` - The mass of the star expressed as a ratio relative to the solar mass.
    /// * `luminosity_in_sols` - The luminosity of the star expressed as a ratio relative to the solar luminosity.
    ///
    /// # Returns
    /// * `age` - The randomly generated current age of the star in years, constrained within realistic limits
    ///   based on the star's expected main sequence lifetime.
    ///
    /// # Example
    /// ```
    /// let star_age = Star::age(1.0, 1.0);
    /// println!("Randomly determined age of a star with solar mass and luminosity: {}", star_age);
    /// ```
    ///
    /// # Notes
    /// - The function first calculates the star's expected main sequence lifetime using `Star::main_seq_life`.
    /// - Depending on this calculated lifetime, the age is generated within different ranges:
    ///   - If the lifetime is at least 6 billion years, the age is randomly chosen between 1 billion and 6 billion years.
    ///   - If the lifetime is between 1 billion and 6 billion years, the age is between 1 billion years and the lifetime.
    ///   - If the lifetime is less than 1 billion years, the age ranges from 1 million to the lifetime, accommodating
    ///     very short-lived stellar phenomena.
    pub fn age(mass_in_sols: f64, luminosity_in_sols: f64) -> f64 {
        let lifetime = Star::main_seq_life(mass_in_sols, luminosity_in_sols);
        let mut rng = rand::thread_rng();

        if lifetime >= 6.0e9 {
            rng.gen_range(1.0e9..=6.0e9)
        } else if lifetime > 1.0e9 {
            rng.gen_range(1.0e9..=lifetime)
        } else {
            rng.gen_range(1.0e6..=lifetime)
        }
    }

    /// Generates a random orbital eccentricity for a celestial object based on its orbital radius.
    /// Eccentricity measures the deviation of an orbit from being circular. An eccentricity of 0 represents
    /// a perfectly circular orbit, whereas values approaching 1 indicate highly elliptical orbits.
    ///
    /// # Parameters
    /// * `orbital_radius_in_au` - The orbital radius of the celestial body in astronomical units (AU).
    ///
    /// # Returns
    /// * `eccentricity` - A double-precision floating-point number representing the orbital eccentricity.
    ///   If the orbital radius is zero (implying the body is at the center or has no discernible orbit),
    ///   the eccentricity is set to zero, representing a circular or undefined orbit.
    ///
    /// # Example
    /// ```
    /// let eccentricity = Star::random_eccentricity(0.5); // For a celestial body at 0.5 AU from the center
    /// println!("Random eccentricity for orbital radius 0.5 AU: {}", eccentricity);
    /// ```
    ///
    /// # Notes
    /// - The function checks if the orbital radius is zero, returning an eccentricity of 0.0 in such cases.
    ///   This condition handles special cases where the orbit might not be defined or is theoretically at the central point.
    /// - For non-zero orbital radii, the function delegates to `AccretionDisk::random_eccentricity`, which should be
    ///   implemented to generate a random eccentricity suitable for typical orbital mechanics in a celestial system.
    pub fn random_eccentricity(orbital_radius_in_au: f64) -> f64 {
        if orbital_radius_in_au == 0.0 {
            0.0
        } else {
            AccretionDisk::random_eccentricity()
        }
    }

    /// Calculates the radius of the ecosphere, also known as the habitable zone, for a star based on its luminosity.
    /// The ecosphere is the region around a star where conditions might be right for liquid water to exist on the surface
    /// of a planet, which is considered crucial for life as we know it. This function uses the square root of the star's
    /// luminosity to estimate the distance at which a planet would need to orbit to potentially support life.
    ///
    /// # Parameters
    /// * `luminosity_in_sols` - The luminosity of the star expressed in solar luminosities (the luminosity of the sun = 1).
    ///
    /// # Returns
    /// * `r_ecosphere` - The radius of the habitable zone in astronomical units (AU), where 1 AU is the average distance
    ///   from the Earth to the Sun.
    ///
    /// # Example
    /// ```
    /// let habitable_zone_radius = Star::r_ecosphere(1.0); // For a star with the luminosity of the sun
    /// println!("Radius of the habitable zone for solar luminosity: {}", habitable_zone_radius);
    /// ```
    ///
    /// # Notes
    /// - This method assumes a simplified model where the square root of the luminosity is proportional to the radius
    ///   of the habitable zone. This is based on the assumption that a star's luminosity affects the amount of radiant
    ///   energy a planet receives, which in turn influences the range of distances at which conditions might support life.
    pub fn r_ecosphere(luminosity_in_sols: f64) -> f64 {
        luminosity_in_sols.sqrt()
    }

    /// Calculates the radius at which a significant greenhouse effect is expected, based on the radius of the ecosphere.
    /// The greenhouse effect radius extends beyond the standard ecosphere radius, factoring in the potential warming
    /// effects due to a planet's atmosphere.
    ///
    /// # Parameters
    /// * `r_ecosphere` - The radius of the ecosphere, typically defined as the distance from a star at which a planet
    ///   can maintain liquid water on its surface under Earth-like conditions, measured in astronomical units (AU).
    ///
    /// # Returns
    /// * `r_greenhouse` - The extended radius accounting for the greenhouse effect, which could potentially allow
    ///   a planet to support liquid water beyond the traditional habitable zone due to atmospheric warming.
    ///
    /// # Example
    /// ```
    /// let r_ecosphere = 1.0; // Ecosphere radius in AU
    /// let r_greenhouse = Star::r_greenhouse(r_ecosphere); // Calculate the greenhouse radius
    /// println!("Greenhouse radius for an ecosphere of 1 AU: {}", r_greenhouse);
    /// ```
    ///
    /// # Notes
    /// - The constant `consts::GREENHOUSE_EFFECT_CONST` is used to scale the ecosphere radius to calculate the greenhouse radius.
    ///   This constant should reflect an empirically or theoretically derived multiplier that considers how much further out
    ///   the capability for liquid water might extend due to atmospheric effects.
    pub fn r_greenhouse(r_ecosphere: f64) -> f64 {
        r_ecosphere * consts::GREENHOUSE_EFFECT_CONST
    }

    pub fn mass_in_sols(max_mass_in_sols: f64) -> f64 {
        random::about(max_mass_in_sols, 0.1)
    }

    pub fn spectral_number(max_mass_in_sols: f64, mass_in_sols: f64, spectral_number: i32) -> i32 {
        let spectral_number_adjustment = (5.0 * (max_mass_in_sols - mass_in_sols) / (max_mass_in_sols)) as i32;
        spectral_number + spectral_number_adjustment
    }

    pub fn calculate_properties(&mut self) {
        self.max_mass_in_sols =
            SpectralInfo::get_max_mass(self.luminosity_class, self.spectral_class, self.spectral_number)
                .expect("Star::from_str failed to find spectral info");
        self.mass_in_sols = Star::mass_in_sols(self.max_mass_in_sols);
        self.spectral_number = Star::spectral_number(self.max_mass_in_sols, self.mass_in_sols, self.spectral_number);
        self.radius_in_au = Star::radius_in_au(self.mass_in_sols, self.luminosity_class, self.spectral_class);
        self.luminosity_in_sols = Star::luminosity_in_sols(self.mass_in_sols, self.luminosity_class);
        self.radius_in_au = Star::radius_in_au(self.mass_in_sols, self.luminosity_class, self.spectral_class);
        self.temperature_in_kelvin = Star::temperature_in_kelvin(self.luminosity_in_sols, self.radius_in_au);
        self.main_seq_life = Star::main_seq_life(self.mass_in_sols, self.luminosity_in_sols);
        self.age = Star::age(self.mass_in_sols, self.luminosity_in_sols);
        self.r_ecosphere = Star::r_ecosphere(self.luminosity_in_sols);
        self.r_greenhouse = Star::r_greenhouse(self.r_ecosphere);
    }

    /// Parses the specification of a star from a command-line input string, typically provided with the "-t" flag.
    ///
    /// This function interprets the input string to extract and construct a `Star` object. The input format
    /// is expected to be "SpectralClassSpectralNumber/LuminosityClassOrbit", where:
    /// - `SpectralClass` is a single character (e.g., 'G').
    /// - `SpectralNumber` is an integer (e.g., 3).
    /// - `LuminosityClass` is a single character (e.g., 'M').
    ///
    /// # Parameters
    /// - `input`: A string slice representing the star's spectral and luminosity information.
    /// - `orbital_radius_in_au`: The orbital radius in astronomical units, not parsed from the input but provided separately.
    ///
    /// # Returns
    /// Returns a `Result
    pub fn from_str(input: &str, orbital_radius_in_au: f64) -> Result<Star, &'static str> {
        // TODO: remove the /orbit part of the input string
        let parts: Vec<&str> = input.split('/').collect();
        if parts.len() != 2 {
            return Err("Input does not match the expected format.");
        }

        let spectral_info = parts[0];
        if spectral_info.len() < 3 {
            return Err("Spectral information is incomplete.");
        }

        let spectral_class_char = spectral_info.chars().next().unwrap();
        let spectral_class =
            SpectralClass::try_from(spectral_class_char).expect("Invalid spectral class {spectral_class_char}");
        let luminosity_id = spectral_info.chars().nth(spectral_info.len() - 1).unwrap();
        let luminosity_class =
            LuminosityClass::try_from(luminosity_id).expect("Invalid luminosity class {luminosity_id}");
        let spectral_number: i32 = match spectral_info[1..spectral_info.len() - 1].parse() {
            Ok(num) => num,
            Err(_) => return Err("Failed to parse spectral number."),
        };

        Ok(Self::new(
            spectral_class,
            luminosity_class,
            spectral_number,
            orbital_radius_in_au,
        ))
    }

    /// Generates a random star (main sequence, dwarf, giant or supergiant), calculating mass, luminosity, etc.
    pub fn random(orbital_radius_in_au: f64) -> Self {
        // According to George Abell's "Exploration of the Universe", (fourth
        // edition), about 90% of all stars in the local neighborhood are main-
        // sequence stars, while about 10% are white dwarfs and less than 1% are
        // giants or supergiants.  This function reflects those percentages.  If
        // you are interested in larger stars, you can always generate them
        //  using the '-t' flag!
        let luminosity_class = match rand::thread_rng().gen_range(1..=100) {
            1..=90 => {
                LuminosityClass::MainSequence // 90% main sequence
            }
            91..=99 => {
                LuminosityClass::WhiteDwarf // 9% white dwarf
            }
            _ => {
                // 1% giants and supergiants
                match rand::thread_rng().gen_range(1..=100) {
                    1..=70 => LuminosityClass::Giant,
                    _ => LuminosityClass::Supergiant,
                }
            }
        };

        // Get a random SpectralInfo entry from one of the STARINFO buckets: main sequence, giant, supergiant, etc.
        let spectral_info = SpectralInfo::get_random(luminosity_class);

        Self::new(
            spectral_info.spec_class,
            luminosity_class,
            spectral_info.spec_num,
            orbital_radius_in_au,
        )
    }

    pub fn new(
        spectral_class: SpectralClass,
        luminosity_class: LuminosityClass,
        spectral_number: i32,
        orbital_radius_in_au: f64,
    ) -> Self {
        let max_mass_in_sols = SpectralInfo::get_max_mass(luminosity_class, spectral_class, spectral_number)
            .expect("Star::from_str failed to find spectral info");
        let mass_in_sols = Star::mass_in_sols(max_mass_in_sols);
        let e = Star::random_eccentricity(orbital_radius_in_au);
        let luminosity_in_sols = Star::luminosity_in_sols(mass_in_sols, luminosity_class);

        let body = Body::new(orbital_radius_in_au, e, mass_in_sols, MassType::Star, 0.0, 0.0);
        let accretion_disk = AccretionDisk::new(body, luminosity_in_sols, orbital_radius_in_au);

        let mut star = Star {
            luminosity_class,
            spectral_class,
            spectral_number,
            temperature_in_kelvin: 0.0,
            orbital_radius_in_au,
            max_mass_in_sols,
            mass_in_sols,
            radius_in_au: 0.0,
            luminosity_in_sols: 0.0,
            age: 0.0,
            main_seq_life: 0.0,
            r_ecosphere: 0.0,
            r_greenhouse: 0.0,
            accretion_disk,
        };

        star.calculate_properties();
        star
    }

    pub fn accrete(&mut self) -> &mut Self {
        self.accretion_disk.accrete();
        self
    }
}
