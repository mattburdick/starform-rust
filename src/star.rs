// src/star.rs

use lazy_static::lazy_static;
use rand::Rng;
use std::collections::HashMap;
use std::fmt;

use crate::accretion_disk::AccretionDisk;
use crate::body::Body;
use crate::types::MassType;
use crate::{consts, random};

lazy_static! {
    #[rustfmt::skip]
        static ref STARINFO: HashMap<LuminosityClass, Vec<SpectralInfo>> = {
        let mut star_info = HashMap::new();

        // Main sequence characteristics
        star_info.insert(LuminosityClass::MainSequence, vec![
            SpectralInfo { spec_class: 'M', spec_num: 9, max_mass: 0.1, percentage: 0.0 },
            SpectralInfo { spec_class: 'M', spec_num: 5, max_mass: 0.2, percentage: 35.0 },
            SpectralInfo { spec_class: 'M', spec_num: 0, max_mass: 0.5, percentage: 36.0 },
            SpectralInfo { spec_class: 'K', spec_num: 5, max_mass: 0.7, percentage: 7.0 },
            SpectralInfo { spec_class: 'K', spec_num: 0, max_mass: 0.8, percentage: 7.0 },
            SpectralInfo { spec_class: 'G', spec_num: 5, max_mass: 0.9, percentage: 3.0 },
            SpectralInfo { spec_class: 'G', spec_num: 0, max_mass: 1.1, percentage: 3.0 },
            SpectralInfo { spec_class: 'F', spec_num: 5, max_mass: 1.3, percentage: 2.0 },
            SpectralInfo { spec_class: 'F', spec_num: 0, max_mass: 1.7, percentage: 1.0 },
            SpectralInfo { spec_class: 'A', spec_num: 5, max_mass: 2.0, percentage: 1.0 },
            SpectralInfo { spec_class: 'A', spec_num: 0, max_mass: 3.2, percentage: 1.0 },
            SpectralInfo { spec_class: 'B', spec_num: 5, max_mass: 6.5, percentage: 1.0 },
            SpectralInfo { spec_class: 'B', spec_num: 0, max_mass: 17.8, percentage: 1.0 },
            SpectralInfo { spec_class: 'O', spec_num: 5, max_mass: 39.8, percentage: 1.0 },
            SpectralInfo { spec_class: 'O', spec_num: 0, max_mass: 60.0, percentage: 1.0 },
        ]);

        // White dwarf characteristics
        star_info.insert(LuminosityClass::WhiteDwarf, vec![
            SpectralInfo { spec_class: 'M', spec_num: 9, max_mass: 0.0, percentage: 0.0 },
            SpectralInfo { spec_class: 'M', spec_num: 5, max_mass: 0.2, percentage: 0.0 },
            SpectralInfo { spec_class: 'M', spec_num: 0, max_mass: 0.4, percentage: 0.0 },
            SpectralInfo { spec_class: 'K', spec_num: 5, max_mass: 0.4, percentage: 1.0 },
            SpectralInfo { spec_class: 'K', spec_num: 0, max_mass: 0.4, percentage: 1.0 },
            SpectralInfo { spec_class: 'G', spec_num: 5, max_mass: 0.5, percentage: 1.0 },
            SpectralInfo { spec_class: 'G', spec_num: 0, max_mass: 0.6, percentage: 1.0 },
            SpectralInfo { spec_class: 'F', spec_num: 5, max_mass: 0.7, percentage: 4.0 },
            SpectralInfo { spec_class: 'F', spec_num: 0, max_mass: 0.8, percentage: 8.0 },
            SpectralInfo { spec_class: 'A', spec_num: 5, max_mass: 1.0, percentage: 28.0 },
            SpectralInfo { spec_class: 'A', spec_num: 0, max_mass: 0.5, percentage: 32.0 },
            SpectralInfo { spec_class: 'B', spec_num: 5, max_mass: 0.4, percentage: 13.0 },
            SpectralInfo { spec_class: 'B', spec_num: 0, max_mass: 0.4, percentage: 9.0 },
            SpectralInfo { spec_class: 'O', spec_num: 5, max_mass: 0.5, percentage: 1.0 },
            SpectralInfo { spec_class: 'O', spec_num: 0, max_mass: 0.7, percentage: 1.0 },
            ]);

        // Giant characteristics
        star_info.insert(LuminosityClass::Giant, vec![
            SpectralInfo { spec_class: 'M', spec_num: 9, max_mass: 8.7, percentage: 0.0 },
            SpectralInfo { spec_class: 'M', spec_num: 5, max_mass: 7.9, percentage: 12.0 },
            SpectralInfo { spec_class: 'M', spec_num: 0, max_mass: 6.3, percentage: 19.0 },
            SpectralInfo { spec_class: 'K', spec_num: 5, max_mass: 5.0, percentage: 26.0 },
            SpectralInfo { spec_class: 'K', spec_num: 0, max_mass: 4.0, percentage: 25.0 },
            SpectralInfo { spec_class: 'G', spec_num: 5, max_mass: 3.2, percentage: 5.0 },
            SpectralInfo { spec_class: 'G', spec_num: 0, max_mass: 2.5, percentage: 4.0 },
            SpectralInfo { spec_class: 'F', spec_num: 5, max_mass: 2.4, percentage: 2.0 },
            SpectralInfo { spec_class: 'F', spec_num: 0, max_mass: 2.5, percentage: 1.0 },
            SpectralInfo { spec_class: 'A', spec_num: 5, max_mass: 2.7, percentage: 1.0 },
            SpectralInfo { spec_class: 'A', spec_num: 0, max_mass: 3.4, percentage: 1.0 },
            SpectralInfo { spec_class: 'B', spec_num: 5, max_mass: 7.0, percentage: 1.0 },
            SpectralInfo { spec_class: 'B', spec_num: 0, max_mass: 30.3, percentage: 1.0 },
            SpectralInfo { spec_class: 'O', spec_num: 5, max_mass: 60.0, percentage: 1.0 },
            SpectralInfo { spec_class: 'O', spec_num: 0, max_mass: 70.0, percentage: 1.0 },
            ]);

        // Supergiant characteristics
        star_info.insert(LuminosityClass::Supergiant, vec![
            SpectralInfo { spec_class: 'M', spec_num: 9, max_mass: 22.3, percentage: 0.0 },
            SpectralInfo { spec_class: 'M', spec_num: 5, max_mass: 19.9, percentage: 12.0 },
            SpectralInfo { spec_class: 'M', spec_num: 0, max_mass: 15.8, percentage: 13.0 },
            SpectralInfo { spec_class: 'K', spec_num: 5, max_mass: 15.0, percentage: 3.0 },
            SpectralInfo { spec_class: 'K', spec_num: 0, max_mass: 12.6, percentage: 4.0 },
            SpectralInfo { spec_class: 'G', spec_num: 5, max_mass: 11.6, percentage: 3.0 },
            SpectralInfo { spec_class: 'G', spec_num: 0, max_mass: 10.0, percentage: 3.0 },
            SpectralInfo { spec_class: 'F', spec_num: 5, max_mass: 11.8, percentage: 8.0 },
            SpectralInfo { spec_class: 'F', spec_num: 0, max_mass: 12.6, percentage: 7.0 },
            SpectralInfo { spec_class: 'A', spec_num: 5, max_mass: 13.2, percentage: 6.0 },
            SpectralInfo { spec_class: 'A', spec_num: 0, max_mass: 15.8, percentage: 6.0 },
            SpectralInfo { spec_class: 'B', spec_num: 5, max_mass: 30.2, percentage: 12.0 },
            SpectralInfo { spec_class: 'B', spec_num: 0, max_mass: 50.1, percentage: 13.0 },
            SpectralInfo { spec_class: 'O', spec_num: 5, max_mass: 70.0, percentage: 4.0 },
            SpectralInfo { spec_class: 'O', spec_num: 0, max_mass: 90.0, percentage: 6.0 },
            ]);

        star_info
    };

}
#[derive(Debug, Copy, Clone)]
struct SpectralInfo {
    spec_class: char,
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
            spec_class: 'G',
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
    fn get_mass(luminosity_class: LuminosityClass, spec_class: char, spec_num: i32) -> Result<f64, &'static str> {
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
    fn to_string(&self) -> &'static str {
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

#[derive(Debug)]
pub struct Star {
    pub luminosity_class: LuminosityClass, // Example: MainSequence, Giant, etc.
    pub spectral_class: char,              // Example: O, B, A, F, G, K, M
    pub spectral_number: i32,              // Example: 3
    pub orbital_radius_in_au: f64,
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

    ///  This is eq. 3.52 from "Astrophysics I" by Bowers and Deeming.
    ///  The mass_in_sols is unitless and is a ratio of the stellar mass to that
    ///  of the Sun.  Both alpha and beta are unitless constants.
    ///  Note that for a main-sequence G3 star like the Sun, this function
    ///  overestimates the luminosity slightly.  It does, however, fit the
    ///  mass-luminosity curve fairly well.
    fn luminosity_in_sols(mass_in_sols: f64, luminosity_class: LuminosityClass) -> Result<f64, &'static str> {
        let luminosity_in_sols: f64;

        let log_mass_ratio = mass_in_sols.log10();
        if log_mass_ratio.is_infinite() {
            return Err("luminosity function: stellar mass ratio leads to infinite logarithm");
        }

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

        Ok(luminosity_in_sols)
    }

    /// This is eq. 3.53 from "Astrophysics I" by Bowers and Deeming.
    /// The mass_ratio is unitless and is a ratio of the stellar mass to that
    /// of the Sun.  The stellar radius returned is in units of AU.
    fn star_radius_in_au(
        mass_in_sols: f64,
        luminosity_class: LuminosityClass,
        spectral_class: char,
    ) -> Result<f64, &'static str> {
        let cool_star = spectral_class == 'K' || spectral_class == 'M';

        let log_mass_ratio = mass_in_sols.log10();
        if log_mass_ratio.is_infinite() {
            return Err("stellar radius function: stellar mass bad");
        }

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

        Ok(radius_in_sols * consts::SOLAR_RADII_PER_AU)
    }

    /// Both the main sequence lifetime and the age returned are in units of
    /// years.  The lifetime passed in is guaranteed to be >= 1 million.
    fn star_age(lifetime: f64) -> f64 {
        let mut rng = rand::thread_rng();

        if lifetime >= 6.0e9 {
            rng.gen_range(1.0e9..=6.0e9)
        } else if lifetime > 1.0e9 {
            rng.gen_range(1.0e9..=lifetime)
        } else {
            rng.gen_range(1.0e6..=lifetime)
        }
    }

    /// Calculate a number of stellar characteristics based on mass and classification
    fn calculate_stellar_stats(
        max_mass: f64,
        luminosity_class: LuminosityClass,
        spectral_class: char,
        spectral_number: i32,
    ) -> (f64, f64, f64, i32, f64, f64, f64, f64) {
        let mass_in_sols = random::about(max_mass, 0.1);
        let luminosity_in_sols =
            Self::luminosity_in_sols(mass_in_sols, luminosity_class).expect("Star::luminosity failed");

        // Adjust the spectral number based on expected stellar mass range for this type of star
        let spectral_number_adjustment = (5.0 * (max_mass - mass_in_sols) / (max_mass)) as i32;
        let spectral_number = spectral_number + spectral_number_adjustment;

        let radius_in_au =
            Self::star_radius_in_au(mass_in_sols, luminosity_class, spectral_class).expect("Star::star_radius failed");

        let main_seq_life = (1.1E10_f64 * (mass_in_sols / luminosity_in_sols)).min(1.0E10_f64);
        let age = Self::star_age(1.0E6_f64.max(1.1E10 * (mass_in_sols / luminosity_in_sols)));
        let r_ecosphere = (luminosity_in_sols).sqrt();
        let r_greenhouse = r_ecosphere * consts::GREENHOUSE_EFFECT_CONST;

        (
            mass_in_sols,
            luminosity_in_sols,
            radius_in_au,
            spectral_number,
            r_ecosphere,
            r_greenhouse,
            age,
            main_seq_life,
        )
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

        let spectral_class = spectral_info.chars().next().unwrap();
        let luminosity_id = spectral_info.chars().nth(spectral_info.len() - 1).unwrap();
        let luminosity_class =
            LuminosityClass::try_from(luminosity_id).expect("Invalid luminosity class {luminosity_id}");
        let spectral_number: i32 = match spectral_info[1..spectral_info.len() - 1].parse() {
            Ok(num) => num,
            Err(_) => return Err("Failed to parse spectral number."),
        };

        let mass_in_sols = SpectralInfo::get_mass(luminosity_class, spectral_class, spectral_number)
            .expect("Star::from_str failed to find spectral info");

        Ok(Self::new(
            spectral_class,
            luminosity_class,
            spectral_number,
            orbital_radius_in_au,
            mass_in_sols,
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
            spectral_info.max_mass,
        )
    }

    pub fn new(
        spectral_class: char,
        luminosity_class: LuminosityClass,
        spectral_number: i32,
        orbital_radius_in_au: f64,
        max_mass_for_class: f64,
    ) -> Self {
        let (mass, luminosity_in_sols, radius_in_au, spectral_number, r_ecosphere, r_greenhouse, age, main_seq_life) =
            Self::calculate_stellar_stats(max_mass_for_class, luminosity_class, spectral_class, spectral_number);
        let e = if orbital_radius_in_au == 0.0 {
            0.0
        } else {
            AccretionDisk::random_eccentricity()
        };

        let star = Body::new(orbital_radius_in_au, e, mass, MassType::Star, 0.0, 0.0);
        let accretion_disk = AccretionDisk::new(star, luminosity_in_sols, orbital_radius_in_au);

        Star {
            spectral_class,
            spectral_number,
            orbital_radius_in_au,
            mass_in_sols: mass,
            age: age,
            luminosity_class,
            main_seq_life,
            r_ecosphere,
            luminosity_in_sols,
            radius_in_au,
            r_greenhouse,
            accretion_disk,
        }
    }

    pub fn accrete(&mut self) -> &mut Self {
        self.accretion_disk.accrete();
        self
    }
}
