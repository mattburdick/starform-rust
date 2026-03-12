use std::fmt;
use std::sync::{Arc, RwLock};

use crate::{
    accretion_disk::AccretionDisk,
    generation::{
        ConfiguredStar, GenerationMode, MultiplicityPreset, SeparationMode, StarCountMode, StarSelectionMode,
        SystemGenerationConfig, MAX_GENERATION_STARS,
    },
    get_log_level, log,
    random::get_random_number,
    star::{LuminosityClass, SpectralClass, Star},
    types::{
        BinaryTopology, HierarchyHostRef, HostedSystem, MassType, PlanetHost, PlanetaryRegion, StableRegion, StarRef,
        StellarHierarchy, StellarHierarchyBinary, StellarHierarchyNode, StellarHierarchyNodeKind, StellarTopology,
    },
};

const MIN_PLANETARY_REGION_WIDTH_AU: f64 = 0.1;
const MIN_HIERARCHICAL_TRIPLE_GAP_RATIO: f64 = 3.0;
use crate::consts::{MIN_PLANETARY_LOG_SPACING_RATIO, MIN_PLANETARY_MUTUAL_HILL_SPACING};
const PLANETARY_MASS_TIE_RELATIVE_TOLERANCE: f64 = 1.0e-6;

#[derive(Debug, Clone, Copy)]
struct HierarchicalOuterOrbit {
    star_index: usize,
    semi_major_axis_au: f64,
    eccentricity: f64,
    gap_ratio: f64,
}

#[derive(Debug, Clone)]
struct PrimaryCenteredHierarchyLayout {
    inner_separation_au: f64,
    inner_eccentricity: f64,
    outer_orbits: Vec<HierarchicalOuterOrbit>,
}

#[derive(Debug, Clone, Copy)]
struct HierarchicalGapFailure {
    star_index: usize,
    inner_apastron_au: f64,
    outer_periastron_au: f64,
    gap_ratio: f64,
}

#[derive(Debug, Clone, Copy)]
struct CircumstellarOuterBounds {
    outer_au: f64,
    stability_outer_au: f64,
    truncation_outer_au: f64,
}

#[derive(Debug, Clone, Copy)]
struct HierarchyAggregate {
    total_mass_solar: f64,
    total_luminosity_solar: f64,
}

#[derive(Debug, Clone)]
struct HierarchyParentBinaryContext {
    parent_node_index: usize,
    sibling_node_index: usize,
    topology: BinaryTopology,
}

#[derive(Debug, Clone)]
pub struct StarSystem {
    pub designation: String,
    pub hosted_systems: Vec<HostedSystem>,
    pub stars: Vec<Star>,
    pub topology: StellarTopology,
    pub planetary_regions: Vec<PlanetaryRegion>,
    pub hierarchy: Option<StellarHierarchy>,
}

// Implement the Display trait for StarSystem
impl fmt::Display for StarSystem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "                         SYSTEM  CHARACTERISTICS")?;
        writeln!(f)?;
        writeln!(f, "Designation: {}", self.display_designation())?;

        for (hosted_system_index, hosted_system) in self.hosted_systems.iter().enumerate() {
            writeln!(f)?;
            if hosted_system_index == 0 {
                writeln!(f, "        PRIMARY HOSTED SYSTEM")?;
            } else {
                writeln!(
                    f,
                    "        HOSTED SYSTEM {}",
                    self.hosted_system_suffix(hosted_system_index)
                )?;
            }
            writeln!(f, "Name: {}", self.hosted_system_designation(hosted_system_index))?;

            for (star_index, star) in hosted_system.stars.iter().enumerate() {
                if star_index == 0 {
                    write!(f, "{}", star)?;
                    if hosted_system.stars.len() > 1 {
                        writeln!(f)?;
                        writeln!(f, "Companion stars present at:")?;
                    }
                } else {
                    writeln!(
                        f,
                        "{:2} {} {:7.3} AU",
                        star_index,
                        star.stellar_classification(),
                        star.a
                    )?;
                }
            }

            match &hosted_system.topology {
                StellarTopology::Single { .. } => {}
                StellarTopology::Binary(binary) => {
                    writeln!(f)?;
                    writeln!(
                        f,
                        "Binary topology: separation {:.3} AU, eccentricity {:.3}",
                        binary.semi_major_axis_au, binary.eccentricity
                    )?;
                }
                StellarTopology::HigherMultiplicity { star_count, .. } => {
                    writeln!(f)?;
                    writeln!(f, "Higher-multiplicity topology: {} stars", star_count)?;
                }
            }
        }

        Ok(())
    }
}

// Implement the StarSystem struct
impl StarSystem {
    pub fn from_stars(stars: Vec<Star>) -> Self {
        let topology = match stars.len() {
            0 | 1 => StellarTopology::Single { primary_index: 0 },
            2 => StellarTopology::Binary(Self::binary_topology_from_stars(&stars)),
            star_count => StellarTopology::HigherMultiplicity {
                primary_index: 0,
                star_count,
            },
        };
        let planetary_regions = Self::collect_hosted_planetary_regions(&stars, &topology, 0);
        let hierarchy = match &topology {
            StellarTopology::Single { .. } => StellarHierarchy::single(0),
            StellarTopology::Binary(binary) => Self::binary_hierarchy(binary),
            StellarTopology::HigherMultiplicity { star_count, .. } => {
                StellarHierarchy::unresolved((0..*star_count).collect())
            }
        };

        Self::from_hosted_systems_with_hierarchy(
            String::new(),
            vec![HostedSystem {
                designation_suffix: None,
                stars,
                topology,
                planetary_regions,
            }],
            hierarchy,
        )
    }

    pub fn from_hosted_systems(designation: impl Into<String>, hosted_systems: Vec<HostedSystem>) -> Self {
        Self::from_hosted_systems_internal(designation, hosted_systems, None)
    }

    pub fn from_hosted_systems_with_hierarchy(
        designation: impl Into<String>,
        hosted_systems: Vec<HostedSystem>,
        hierarchy: StellarHierarchy,
    ) -> Self {
        Self::from_hosted_systems_internal(designation, hosted_systems, Some(hierarchy))
    }

    fn from_hosted_systems_internal(
        designation: impl Into<String>,
        hosted_systems: Vec<HostedSystem>,
        hierarchy: Option<StellarHierarchy>,
    ) -> Self {
        let designation = designation.into();
        let stars = hosted_systems
            .iter()
            .flat_map(|hosted_system| hosted_system.stars.iter().cloned())
            .collect::<Vec<_>>();
        let planetary_regions = hosted_systems
            .iter()
            .flat_map(|hosted_system| hosted_system.planetary_regions.iter().cloned())
            .collect::<Vec<_>>();
        let topology = Self::legacy_root_topology(&hosted_systems, stars.len());

        StarSystem {
            designation,
            hosted_systems,
            topology,
            planetary_regions,
            stars,
            hierarchy,
        }
    }

    pub fn with_designation(mut self, designation: impl Into<String>) -> Self {
        self.designation = designation.into();
        self
    }

    pub fn ensure_designation(&mut self, designation: impl Into<String>) {
        if self.designation.trim().is_empty() {
            self.designation = designation.into();
        }
    }

    pub fn display_designation(&self) -> String {
        if self.designation.trim().is_empty() {
            "Unnamed".to_string()
        } else {
            self.designation.clone()
        }
    }

    pub fn hosted_system(&self, hosted_system_index: usize) -> Option<&HostedSystem> {
        self.hosted_systems.get(hosted_system_index)
    }

    pub fn star(&self, star_ref: StarRef) -> Option<&Star> {
        self.hosted_system(star_ref.hosted_system_index)
            .and_then(|hosted_system| hosted_system.stars.get(star_ref.star_index))
    }

    pub fn primary_star_ref(&self) -> Option<StarRef> {
        self.hosted_system(0).and_then(|hosted_system| {
            hosted_system.stars.first().map(|_| StarRef {
                hosted_system_index: 0,
                star_index: 0,
            })
        })
    }

    pub fn hosted_system_designation(&self, hosted_system_index: usize) -> String {
        if hosted_system_index == 0 {
            return self.display_designation();
        }

        format!(
            "{} {}",
            self.display_designation(),
            self.hosted_system_suffix(hosted_system_index)
        )
    }

    pub fn hosted_system_suffix(&self, hosted_system_index: usize) -> String {
        if hosted_system_index == 0 {
            return self
                .hosted_system(0)
                .and_then(|hosted_system| hosted_system.designation_suffix.clone())
                .unwrap_or_default();
        }

        self.hosted_system(hosted_system_index)
            .and_then(|hosted_system| hosted_system.designation_suffix.clone())
            .unwrap_or_else(|| Self::default_hosted_system_suffix(hosted_system_index))
    }

    fn legacy_root_topology(hosted_systems: &[HostedSystem], total_star_count: usize) -> StellarTopology {
        if let Some(hosted_system) = hosted_systems.first() {
            if hosted_systems.len() == 1 {
                return hosted_system.topology.clone();
            }
        }

        StellarTopology::HigherMultiplicity {
            primary_index: 0,
            star_count: total_star_count,
        }
    }

    fn default_hosted_system_suffix(hosted_system_index: usize) -> String {
        let mut ordinal = hosted_system_index;
        let mut suffix = String::new();

        loop {
            let letter = ((ordinal % 26) as u8 + b'A') as char;
            suffix.insert(0, letter);
            if ordinal < 26 {
                break;
            }
            ordinal = (ordinal / 26) - 1;
        }

        suffix
    }

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
    /// ```rust
    /// use starform_rust::star_system::StarSystem;
    ///
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

            match star_count {
                1 => {
                    let star = Star::random(0.0);
                    log!(
                        log_level,
                        1,
                        "Random star: {:.2} solar masses - {}",
                        star.mass_in_sols,
                        star.stellar_classification()
                    );
                    Self::generate_single_star_system(star, GenerationMode::Aggregation)
                }
                2 => Self::generate_binary_star_system(log_level, GenerationMode::Aggregation),
                3 => Self::generate_random_triple_star_system(log_level, GenerationMode::Aggregation),
                _ => Self::generate_random_higher_multiplicity_star_system(
                    log_level,
                    star_count,
                    GenerationMode::Aggregation,
                ),
            }
        } else {
            // Create a star system with one star based on the description in the "-t" flag
            let orbital_radius_in_au = 0.0;
            let star = match Star::from_str(star_type, orbital_radius_in_au) {
                Ok(value) => value,
                Err(err) => panic!("Error: {}", err),
            };

            log!(
                log_level,
                1,
                "Star: {:.2} solar masses - {}",
                star.mass_in_sols,
                star.stellar_classification()
            );
            Self::generate_single_star_system(star, GenerationMode::Aggregation)
        }
    }

    /// Create configured stars from a generation config **without** running
    /// accretion.  Each returned star has an initialised (but un-accreted)
    /// `AccretionDisk`.  The RNG is advanced exactly as `new_with_config` would
    /// up to (but not including) the accretion step, so callers that later run
    /// accretion on each star and re-assemble via `from_stars` will get the
    /// same result as the one-shot path.
    ///
    /// This is the entry-point used by the animated-generation UI: it creates
    /// `AccretionStepper`s from the returned stars' disks and drives them one
    /// step at a time.
    pub fn prepare_stars_from_config(generation_config: &SystemGenerationConfig) -> Vec<Star> {
        let log_level = *get_log_level!();
        let stars = Self::build_stars_from_config(generation_config);
        Self::log_generated_stars(log_level, &stars);
        Self::log_multi_star_layout(log_level, &stars);
        stars
    }

    fn build_stars_from_config(generation_config: &SystemGenerationConfig) -> Vec<Star> {
        let star_count = Self::resolve_star_count(generation_config);
        let companion_separations_au = Self::resolve_companion_separations_au(generation_config, star_count);

        let mut stars = Vec::with_capacity(star_count);
        for star_index in 0..star_count {
            let orbital_radius_in_au = if star_index == 0 {
                0.0
            } else {
                companion_separations_au[star_index - 1]
            };
            let mut star = Self::build_configured_star(&generation_config.stars[star_index], orbital_radius_in_au);
            let eccentricity = Star::random_eccentricity(orbital_radius_in_au);
            Self::configure_star_orbit(&mut star, orbital_radius_in_au, eccentricity);
            stars.push(star);
        }
        stars
    }

    pub fn new_with_config(generation_config: &SystemGenerationConfig) -> Self {
        let log_level = *get_log_level!();
        let generation_mode = generation_config.generation_mode;
        let mut stars = Self::build_stars_from_config(generation_config);

        Self::log_generated_stars(log_level, &stars);
        Self::log_multi_star_layout(log_level, &stars);

        match stars.len() {
            1 => Self::generate_single_star_system(stars.remove(0), generation_mode),
            2 => {
                let binary = BinaryTopology {
                    primary_index: 0,
                    secondary_index: 1,
                    semi_major_axis_au: stars[1].a,
                    eccentricity: stars[1].body.e,
                    total_mass_solar: stars[0].mass_in_sols + stars[1].mass_in_sols,
                    mass_ratio: if (stars[0].mass_in_sols + stars[1].mass_in_sols) > 0.0 {
                        stars[1].mass_in_sols / (stars[0].mass_in_sols + stars[1].mass_in_sols)
                    } else {
                        0.0
                    },
                };
                let topology = StellarTopology::Binary(binary.clone());
                let planetary_regions = Self::build_binary_planetary_regions(&mut stars, &topology, 0, generation_mode);
                Self::from_hosted_systems_with_hierarchy(
                    String::new(),
                    vec![HostedSystem {
                        designation_suffix: None,
                        stars,
                        topology,
                        planetary_regions,
                    }],
                    Self::binary_hierarchy(&binary),
                )
            }
            3 => Self::generate_triple_star_system(log_level, stars, generation_mode),
            _ => Self::generate_higher_multiplicity_star_system(log_level, stars, generation_mode),
        }
    }

    fn generate_random_triple_star_system(log_level: u8, generation_mode: GenerationMode) -> Self {
        let inner_separation_au = Self::sample_observed_binary_separation_au();
        let outer_separation_au = Self::sample_outer_companion_separation_au(inner_separation_au);
        let orbital_radii_au = [0.0, inner_separation_au, outer_separation_au];
        let mut stars = orbital_radii_au
            .into_iter()
            .map(|orbital_radius_in_au| {
                let mut star = Star::random(orbital_radius_in_au);
                let eccentricity = Star::random_eccentricity(orbital_radius_in_au);
                Self::configure_star_orbit(&mut star, orbital_radius_in_au, eccentricity);
                star
            })
            .collect::<Vec<_>>();

        Self::log_generated_stars(log_level, &stars);
        Self::log_multi_star_layout(log_level, &stars);
        Self::generate_triple_star_system(log_level, std::mem::take(&mut stars), generation_mode)
    }

    fn generate_single_star_system(mut star: Star, generation_mode: GenerationMode) -> Self {
        star.accrete_with_mode(generation_mode);
        let stars = vec![star];
        let topology = StellarTopology::Single { primary_index: 0 };
        let planetary_regions = Self::collect_hosted_planetary_regions(&stars, &topology, 0);

        Self::from_hosted_systems_with_hierarchy(
            String::new(),
            vec![HostedSystem {
                designation_suffix: None,
                stars,
                topology,
                planetary_regions,
            }],
            StellarHierarchy::single(0),
        )
    }

    fn generate_random_higher_multiplicity_star_system(
        log_level: u8,
        star_count: usize,
        generation_mode: GenerationMode,
    ) -> Self {
        let mut companion_separations_au = Vec::with_capacity(star_count.saturating_sub(1));
        if star_count > 1 {
            let first_separation_au = Self::sample_observed_binary_separation_au();
            companion_separations_au.push(first_separation_au);
            while companion_separations_au.len() < star_count - 1 {
                let previous = *companion_separations_au.last().unwrap_or(&first_separation_au);
                companion_separations_au.push(Self::sample_outer_companion_separation_au(previous));
            }
        }

        let mut stars: Vec<Star> = Vec::with_capacity(star_count);
        for index in 0..star_count {
            let orbital_radius_in_au = if index == 0 {
                0.0
            } else {
                companion_separations_au[index - 1]
            };

            let mut star = Star::random(orbital_radius_in_au);
            let eccentricity = Star::random_eccentricity(orbital_radius_in_au);
            Self::configure_star_orbit(&mut star, orbital_radius_in_au, eccentricity);
            log!(
                log_level,
                1,
                "Random star: {:.2} solar masses - {}",
                star.mass_in_sols,
                star.stellar_classification()
            );
            stars.push(star);
        }

        Self::log_multi_star_layout(log_level, &stars);
        Self::generate_higher_multiplicity_star_system(log_level, stars, generation_mode)
    }

    fn generate_binary_star_system(log_level: u8, generation_mode: GenerationMode) -> Self {
        let mut stars = vec![Star::random(0.0), Star::random(0.0)];
        for star in &stars {
            log!(
                log_level,
                1,
                "Random star: {:.2} solar masses - {}",
                star.mass_in_sols,
                star.stellar_classification()
            );
        }

        stars.sort_by(|left, right| {
            right
                .mass_in_sols
                .partial_cmp(&left.mass_in_sols)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let binary_separation_au = get_random_number(1.0..=150.0);
        let binary_eccentricity = Star::random_eccentricity(binary_separation_au);
        Self::configure_star_orbit(&mut stars[0], 0.0, 0.0);
        Self::configure_star_orbit(&mut stars[1], binary_separation_au, binary_eccentricity);
        Self::log_multi_star_layout(log_level, &stars);

        let binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: binary_separation_au,
            eccentricity: binary_eccentricity,
            total_mass_solar: stars[0].mass_in_sols + stars[1].mass_in_sols,
            mass_ratio: stars[1].mass_in_sols / (stars[0].mass_in_sols + stars[1].mass_in_sols),
        };
        let topology = StellarTopology::Binary(binary.clone());
        let planetary_regions = Self::build_binary_planetary_regions(&mut stars, &topology, 0, generation_mode);

        Self::from_hosted_systems_with_hierarchy(
            String::new(),
            vec![HostedSystem {
                designation_suffix: None,
                stars,
                topology,
                planetary_regions,
            }],
            Self::binary_hierarchy(&binary),
        )
    }

    fn generate_triple_star_system(log_level: u8, mut stars: Vec<Star>, generation_mode: GenerationMode) -> Self {
        let layout = match Self::resolve_primary_centered_hierarchy(&stars) {
            Ok(layout) if layout.outer_orbits.len() == 1 => layout,
            Ok(_) => {
                return Self::build_unresolved_multi_star_system(String::new(), std::mem::take(&mut stars));
            }
            Err(failure) => {
                log!(
                    log_level,
                    1,
                    "DEBUG-MULTISTAR rejecting hierarchical triple planet generation: current sampled primary-centered layout is not clearly hierarchical (inner_apastron={:.3} AU, outer_periastron={:.3} AU, gap_ratio={:.3}, minimum={:.3}); returning stars without planetary regions to avoid misleading overlapping disks",
                    failure.inner_apastron_au,
                    failure.outer_periastron_au,
                    failure.gap_ratio,
                    MIN_HIERARCHICAL_TRIPLE_GAP_RATIO
                );
                return Self::build_unresolved_multi_star_system(String::new(), std::mem::take(&mut stars));
            }
        };
        let outer_orbit = layout.outer_orbits[0];

        log!(
            log_level,
            1,
            "DEBUG-MULTISTAR resolved hierarchical triple inner_pair=(0,1) a={:.3} AU e={:.3} outer_tertiary=2 a={:.3} AU e={:.3} gap_ratio={:.3}",
            layout.inner_separation_au,
            layout.inner_eccentricity,
            outer_orbit.semi_major_axis_au,
            outer_orbit.eccentricity,
            outer_orbit.gap_ratio
        );

        Self::build_primary_centered_hierarchical_system(String::new(), stars, &layout, log_level, generation_mode)
    }

    fn generate_higher_multiplicity_star_system(
        log_level: u8,
        stars: Vec<Star>,
        generation_mode: GenerationMode,
    ) -> Self {
        let star_count = stars.len();
        let layout = match Self::resolve_primary_centered_hierarchy(&stars) {
            Ok(layout) => layout,
            Err(failure) => {
                log!(
                    log_level,
                    1,
                    "DEBUG-MULTISTAR rejecting higher-multiplicity hierarchy resolution for {} stars: companion {} does not maintain a clear hierarchical gap (inner_apastron={:.3} AU, outer_periastron={:.3} AU, gap_ratio={:.3}, minimum={:.3}); returning stars without planetary regions",
                    star_count,
                    failure.star_index,
                    failure.inner_apastron_au,
                    failure.outer_periastron_au,
                    failure.gap_ratio,
                    MIN_HIERARCHICAL_TRIPLE_GAP_RATIO
                );
                return Self::build_unresolved_multi_star_system(String::new(), stars);
            }
        };

        let outer_summary = layout
            .outer_orbits
            .iter()
            .map(|orbit| {
                format!(
                    "{}:{:.3} AU e={:.3} gap_ratio={:.3}",
                    orbit.star_index, orbit.semi_major_axis_au, orbit.eccentricity, orbit.gap_ratio
                )
            })
            .collect::<Vec<_>>();
        log!(
            log_level,
            1,
            "DEBUG-MULTISTAR resolved primary-centered {}-star hierarchy inner_pair=(0,1) a={:.3} AU e={:.3} outers=[{}]; Slice 3 will derive conservative regions only for directly justified star and binary hosts",
            star_count,
            layout.inner_separation_au,
            layout.inner_eccentricity,
            outer_summary.join(", ")
        );

        Self::build_primary_centered_hierarchical_system(String::new(), stars, &layout, log_level, generation_mode)
    }

    fn resolve_star_count(generation_config: &SystemGenerationConfig) -> usize {
        match generation_config.star_count_mode {
            StarCountMode::Exact => generation_config.star_count.clamp(1, MAX_GENERATION_STARS),
            StarCountMode::Sampled => Self::sample_star_count(generation_config.multiplicity_preset),
        }
    }

    fn sample_star_count(multiplicity_preset: MultiplicityPreset) -> usize {
        match multiplicity_preset {
            MultiplicityPreset::Custom | MultiplicityPreset::Observed => match get_random_number(1..=100) {
                1..=54 => 1,
                55..=87 => 2,
                88..=96 => 3,
                _ => 4,
            },
            MultiplicityPreset::SingleHeavy => match get_random_number(1..=100) {
                1..=74 => 1,
                75..=94 => 2,
                95..=99 => 3,
                _ => 4,
            },
            MultiplicityPreset::BinaryHeavy => match get_random_number(1..=100) {
                1..=28 => 1,
                29..=78 => 2,
                79..=94 => 3,
                _ => 4,
            },
            MultiplicityPreset::HigherMultiplicity => match get_random_number(1..=100) {
                1..=12 => 1,
                13..=40 => 2,
                41..=74 => 3,
                _ => 4,
            },
        }
    }

    fn resolve_companion_separations_au(generation_config: &SystemGenerationConfig, star_count: usize) -> Vec<f64> {
        if star_count <= 1 {
            return Vec::new();
        }

        let mut separations_au = Vec::with_capacity(star_count - 1);
        let first_separation_au = match generation_config.separation_mode {
            SeparationMode::Sampled => Self::sample_observed_binary_separation_au(),
            SeparationMode::Exact => generation_config.binary_separation_au.clamp(0.05, 5_000.0),
        };
        separations_au.push(first_separation_au);

        while separations_au.len() < star_count - 1 {
            let previous = *separations_au.last().unwrap_or(&first_separation_au);
            let next = match generation_config.separation_mode {
                SeparationMode::Sampled => Self::sample_outer_companion_separation_au(previous),
                SeparationMode::Exact => (previous * 4.0).clamp(0.2, 5_000.0),
            };
            separations_au.push(next);
        }

        separations_au
    }

    fn sample_observed_binary_separation_au() -> f64 {
        match get_random_number(1..=1000) {
            1..=80 => get_random_number(0.02..=0.2),
            81..=320 => get_random_number(0.2..=5.0),
            321..=670 => get_random_number(5.0..=60.0),
            671..=900 => get_random_number(60.0..=300.0),
            _ => get_random_number(300.0..=2_000.0),
        }
    }

    fn sample_outer_companion_separation_au(previous_separation_au: f64) -> f64 {
        let max_separation_au = 5_000.0;
        let min_separation_au = (previous_separation_au * 1.8).min(max_separation_au);
        if min_separation_au >= max_separation_au {
            return max_separation_au;
        }

        let multiplier = match get_random_number(1..=100) {
            1..=25 => get_random_number(2.5..=4.0),
            26..=75 => get_random_number(4.0..=8.0),
            _ => get_random_number(8.0..=16.0),
        };
        (previous_separation_au * multiplier).clamp(min_separation_au, max_separation_au)
    }

    fn build_configured_star(configured_star: &ConfiguredStar, orbital_radius_in_au: f64) -> Star {
        match configured_star.selection_mode {
            StarSelectionMode::Explicit => Star::new(
                configured_star.spectral_class,
                configured_star.luminosity_class,
                configured_star.spectral_number.clamp(0, 9),
                orbital_radius_in_au,
            ),
            StarSelectionMode::RandomObserved => Self::sample_random_observed_star(orbital_radius_in_au),
        }
    }

    fn sample_random_observed_star(orbital_radius_in_au: f64) -> Star {
        let luminosity_class = Self::sample_observed_luminosity_class();
        let spectral_class = Self::sample_observed_spectral_class(luminosity_class);
        let spectral_number = get_random_number(0..=9);
        Star::new(spectral_class, luminosity_class, spectral_number, orbital_radius_in_au)
    }

    fn sample_observed_luminosity_class() -> LuminosityClass {
        match get_random_number(1..=10_000) {
            1..=9_250 => LuminosityClass::MainSequence,
            9_251..=9_900 => LuminosityClass::WhiteDwarf,
            9_901..=9_960 => LuminosityClass::Subgiant,
            9_961..=9_990 => LuminosityClass::Giant,
            9_991..=9_997 => LuminosityClass::BrightGiant,
            9_998..=9_999 => LuminosityClass::Supergiant,
            _ => LuminosityClass::Hypergiant,
        }
    }

    fn sample_observed_spectral_class(luminosity_class: LuminosityClass) -> SpectralClass {
        match luminosity_class {
            LuminosityClass::MainSequence => match get_random_number(1..=1000) {
                1..=765 => SpectralClass::M,
                766..=885 => SpectralClass::K,
                886..=955 => SpectralClass::G,
                956..=985 => SpectralClass::F,
                986..=996 => SpectralClass::A,
                997..=999 => SpectralClass::B,
                _ => SpectralClass::O,
            },
            LuminosityClass::WhiteDwarf => match get_random_number(1..=100) {
                1..=58 => SpectralClass::A,
                59..=78 => SpectralClass::B,
                79..=88 => SpectralClass::F,
                89..=94 => SpectralClass::G,
                95..=98 => SpectralClass::K,
                _ => SpectralClass::O,
            },
            LuminosityClass::Subgiant => match get_random_number(1..=100) {
                1..=24 => SpectralClass::K,
                25..=48 => SpectralClass::G,
                49..=68 => SpectralClass::F,
                69..=83 => SpectralClass::A,
                84..=92 => SpectralClass::M,
                93..=98 => SpectralClass::B,
                _ => SpectralClass::O,
            },
            LuminosityClass::Giant => match get_random_number(1..=100) {
                1..=44 => SpectralClass::K,
                45..=79 => SpectralClass::M,
                80..=89 => SpectralClass::G,
                90..=94 => SpectralClass::F,
                95..=97 => SpectralClass::A,
                98..=99 => SpectralClass::B,
                _ => SpectralClass::O,
            },
            LuminosityClass::BrightGiant => match get_random_number(1..=100) {
                1..=10 => SpectralClass::M,
                11..=22 => SpectralClass::K,
                23..=38 => SpectralClass::G,
                39..=58 => SpectralClass::F,
                59..=75 => SpectralClass::A,
                76..=92 => SpectralClass::B,
                _ => SpectralClass::O,
            },
            LuminosityClass::Supergiant | LuminosityClass::Hypergiant => match get_random_number(1..=100) {
                1..=10 => SpectralClass::M,
                11..=18 => SpectralClass::K,
                19..=28 => SpectralClass::G,
                29..=44 => SpectralClass::F,
                45..=60 => SpectralClass::A,
                61..=84 => SpectralClass::B,
                _ => SpectralClass::O,
            },
        }
    }

    fn log_generated_stars(log_level: u8, stars: &[Star]) {
        for star in stars {
            log!(
                log_level,
                1,
                "Random/configured star: {:.2} solar masses - {}",
                star.mass_in_sols,
                star.stellar_classification()
            );
        }
    }

    fn log_multi_star_layout(log_level: u8, stars: &[Star]) {
        if stars.len() <= 1 {
            return;
        }

        let companions = stars
            .iter()
            .enumerate()
            .skip(1)
            .map(|(star_index, star)| {
                format!(
                    "{}:{:.3} AU ({}, e={:.3})",
                    star_index,
                    star.a,
                    star.stellar_classification(),
                    star.body.e
                )
            })
            .collect::<Vec<_>>();

        log!(
            log_level,
            1,
            "DEBUG-MULTISTAR companion separations relative_to_primary=[{}]",
            companions.join(", ")
        );

        if stars.len() >= 4 {
            log!(
                log_level,
                1,
                "DEBUG-MULTISTAR higher-multiplicity generation will attempt topology-first hierarchy resolution for 4+ stars before any planet generation"
            );
        }
    }

    fn collect_hosted_planetary_regions(
        stars: &[Star],
        topology: &StellarTopology,
        hosted_system_index: usize,
    ) -> Vec<PlanetaryRegion> {
        stars
            .iter()
            .enumerate()
            .filter_map(|(star_index, star)| {
                let accretion_disk = star.body.accretion_disk.as_ref()?.clone();
                let disk = match accretion_disk.read() {
                    Ok(guard) => guard,
                    Err(poisoned) => poisoned.into_inner(),
                };
                Some(Self::build_planetary_region(
                    StableRegion {
                        host: PlanetHost::Circumstellar {
                            star_ref: StarRef {
                                hosted_system_index,
                                star_index,
                            },
                        },
                        hierarchy_host: Self::simple_topology_star_hierarchy_host(topology, star_index),
                        inner_au: disk.planet_inner_bound,
                        outer_au: disk.planet_outer_bound,
                        stability_inner_au: disk.planet_inner_bound,
                        stability_outer_au: disk.planet_outer_bound,
                        truncation_inner_au: disk.planet_inner_bound,
                        truncation_outer_au: disk.planet_outer_bound,
                    },
                    star.mass_in_sols,
                    star.luminosity_in_sols,
                    disk.bodies.clone(),
                ))
            })
            .collect()
    }

    fn simple_topology_star_hierarchy_host(topology: &StellarTopology, star_index: usize) -> Option<HierarchyHostRef> {
        match topology {
            StellarTopology::Single { .. } if star_index == 0 => Some(HierarchyHostRef::Star {
                node_index: 0,
                global_star_index: 0,
            }),
            StellarTopology::Binary(_) if star_index < 2 => Some(HierarchyHostRef::Star {
                node_index: star_index,
                global_star_index: star_index,
            }),
            _ => None,
        }
    }

    fn hierarchy_star_host(node_index: usize, global_star_index: usize) -> HierarchyHostRef {
        HierarchyHostRef::Star {
            node_index,
            global_star_index,
        }
    }

    fn hierarchy_binary_host(node_index: usize, binary: &StellarHierarchyBinary) -> HierarchyHostRef {
        HierarchyHostRef::Binary {
            node_index,
            primary_child_index: binary.primary_child_index,
            secondary_child_index: binary.secondary_child_index,
        }
    }

    fn build_primary_centered_hierarchical_system(
        designation: impl Into<String>,
        mut stars: Vec<Star>,
        layout: &PrimaryCenteredHierarchyLayout,
        log_level: u8,
        generation_mode: GenerationMode,
    ) -> Self {
        let hierarchy = Self::primary_centered_hierarchy(layout);

        if let Some(primary) = stars.get_mut(0) {
            Self::configure_star_orbit(primary, 0.0, 0.0);
        }
        if let Some(secondary) = stars.get_mut(1) {
            Self::configure_star_orbit(secondary, layout.inner_separation_au, layout.inner_eccentricity);
        }
        for outer_orbit in &layout.outer_orbits {
            if let Some(outer_star) = stars.get_mut(outer_orbit.star_index) {
                Self::configure_star_orbit(outer_star, outer_orbit.semi_major_axis_au, outer_orbit.eccentricity);
            }
        }

        let hosted_system_count = layout.outer_orbits.len() + 1;
        let mut planetary_regions_by_hosted_system = vec![Vec::new(); hosted_system_count];
        for region in Self::build_hierarchy_planetary_regions(log_level, &mut stars, &hierarchy, generation_mode) {
            let hosted_system_index = region.stable_region.host.hosted_system_index();
            if let Some(hosted_system_regions) = planetary_regions_by_hosted_system.get_mut(hosted_system_index) {
                hosted_system_regions.push(region);
            }
        }
        for hosted_system_regions in &mut planetary_regions_by_hosted_system {
            Self::sort_planetary_regions(hosted_system_regions);
        }

        let inner_primary_mass = stars.first().map(|star| star.mass_in_sols).unwrap_or(0.0);
        let inner_secondary_mass = stars.get(1).map(|star| star.mass_in_sols).unwrap_or(0.0);
        let inner_total_mass_solar = inner_primary_mass + inner_secondary_mass;
        let inner_binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: layout.inner_separation_au,
            eccentricity: layout.inner_eccentricity,
            total_mass_solar: inner_total_mass_solar,
            mass_ratio: if inner_total_mass_solar > 0.0 {
                inner_secondary_mass / inner_total_mass_solar
            } else {
                0.0
            },
        };

        let inner_primary = stars[0].clone();
        let inner_secondary = stars[1].clone();

        let mut hosted_systems = vec![HostedSystem {
            designation_suffix: None,
            stars: vec![inner_primary, inner_secondary],
            topology: StellarTopology::Binary(inner_binary),
            planetary_regions: std::mem::take(&mut planetary_regions_by_hosted_system[0]),
        }];

        for outer_orbit in &layout.outer_orbits {
            let hosted_system_index = hosted_systems.len();
            let outer_star = stars[outer_orbit.star_index].clone();
            hosted_systems.push(HostedSystem {
                designation_suffix: Some(Self::default_hosted_system_suffix(hosted_system_index)),
                stars: vec![outer_star],
                topology: StellarTopology::Single { primary_index: 0 },
                planetary_regions: std::mem::take(&mut planetary_regions_by_hosted_system[hosted_system_index]),
            });
        }

        Self::from_hosted_systems_with_hierarchy(designation, hosted_systems, hierarchy)
    }

    fn build_unresolved_multi_star_system(designation: impl Into<String>, mut stars: Vec<Star>) -> Self {
        let star_count = stars.len();
        for star in &mut stars {
            let empty_disk = Self::empty_unbounded_circumstellar_disk(star);
            star.body.accretion_disk = Some(Arc::new(RwLock::new(empty_disk)));
        }

        Self::from_hosted_systems_with_hierarchy(
            designation,
            vec![HostedSystem {
                designation_suffix: None,
                stars,
                topology: StellarTopology::HigherMultiplicity {
                    primary_index: 0,
                    star_count,
                },
                planetary_regions: Vec::new(),
            }],
            StellarHierarchy::unresolved((0..star_count).collect()),
        )
    }

    fn resolve_primary_centered_hierarchy(
        stars: &[Star],
    ) -> Result<PrimaryCenteredHierarchyLayout, HierarchicalGapFailure> {
        let inner = stars.get(1).ok_or(HierarchicalGapFailure {
            star_index: 1,
            inner_apastron_au: 0.0,
            outer_periastron_au: 0.0,
            gap_ratio: 0.0,
        })?;
        if inner.a <= 0.0 {
            return Err(HierarchicalGapFailure {
                star_index: 1,
                inner_apastron_au: 0.0,
                outer_periastron_au: 0.0,
                gap_ratio: 0.0,
            });
        }

        let mut current_apastron_au = inner.a * (1.0 + inner.body.e.max(0.0));
        if current_apastron_au <= 0.0 {
            return Err(HierarchicalGapFailure {
                star_index: 1,
                inner_apastron_au: current_apastron_au,
                outer_periastron_au: 0.0,
                gap_ratio: 0.0,
            });
        }

        let mut previous_separation_au = inner.a;
        let mut outer_orbits = Vec::with_capacity(stars.len().saturating_sub(2));
        for (star_index, outer) in stars.iter().enumerate().skip(2) {
            let outer_periastron_au = outer.a * (1.0 - outer.body.e.clamp(0.0, 0.999));
            let gap_ratio = if current_apastron_au > 0.0 {
                outer_periastron_au / current_apastron_au
            } else {
                0.0
            };
            if outer.a <= previous_separation_au
                || outer_periastron_au <= 0.0
                || gap_ratio < MIN_HIERARCHICAL_TRIPLE_GAP_RATIO
            {
                return Err(HierarchicalGapFailure {
                    star_index,
                    inner_apastron_au: current_apastron_au,
                    outer_periastron_au,
                    gap_ratio,
                });
            }

            outer_orbits.push(HierarchicalOuterOrbit {
                star_index,
                semi_major_axis_au: outer.a,
                eccentricity: outer.body.e,
                gap_ratio,
            });
            current_apastron_au = outer.a * (1.0 + outer.body.e.max(0.0));
            previous_separation_au = outer.a;
        }

        Ok(PrimaryCenteredHierarchyLayout {
            inner_separation_au: inner.a,
            inner_eccentricity: inner.body.e,
            outer_orbits,
        })
    }

    fn binary_topology_from_stars(stars: &[Star]) -> BinaryTopology {
        let primary_mass = stars.first().map(|star| star.mass_in_sols).unwrap_or(0.0);
        let secondary_mass = stars.get(1).map(|star| star.mass_in_sols).unwrap_or(0.0);
        let total_mass = primary_mass + secondary_mass;
        let semi_major_axis_au = stars.get(1).map(|star| star.a).unwrap_or(0.0);
        let eccentricity = stars.get(1).map(|star| star.body.e).unwrap_or(0.0);

        BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au,
            eccentricity,
            total_mass_solar: total_mass,
            mass_ratio: if total_mass > 0.0 {
                secondary_mass / total_mass
            } else {
                0.0
            },
        }
    }

    fn binary_hierarchy(binary: &BinaryTopology) -> StellarHierarchy {
        let root_index = 2;
        StellarHierarchy {
            root_index,
            nodes: vec![
                StellarHierarchyNode {
                    parent_index: Some(root_index),
                    kind: StellarHierarchyNodeKind::Star {
                        star_index: binary.primary_index,
                    },
                },
                StellarHierarchyNode {
                    parent_index: Some(root_index),
                    kind: StellarHierarchyNodeKind::Star {
                        star_index: binary.secondary_index,
                    },
                },
                StellarHierarchyNode {
                    parent_index: None,
                    kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                        primary_child_index: 0,
                        secondary_child_index: 1,
                        semi_major_axis_au: binary.semi_major_axis_au,
                        eccentricity: binary.eccentricity,
                    }),
                },
            ],
        }
    }

    fn primary_centered_hierarchy(layout: &PrimaryCenteredHierarchyLayout) -> StellarHierarchy {
        let initial_root_index = 2;
        let mut nodes = vec![
            StellarHierarchyNode {
                parent_index: Some(initial_root_index),
                kind: StellarHierarchyNodeKind::Star { star_index: 0 },
            },
            StellarHierarchyNode {
                parent_index: Some(initial_root_index),
                kind: StellarHierarchyNodeKind::Star { star_index: 1 },
            },
            StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 0,
                    secondary_child_index: 1,
                    semi_major_axis_au: layout.inner_separation_au,
                    eccentricity: layout.inner_eccentricity,
                }),
            },
        ];
        let mut current_root_index = initial_root_index;

        for outer_orbit in &layout.outer_orbits {
            let star_node_index = nodes.len();
            nodes.push(StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Star {
                    star_index: outer_orbit.star_index,
                },
            });

            let new_root_index = nodes.len();
            nodes[current_root_index].parent_index = Some(new_root_index);
            nodes[star_node_index].parent_index = Some(new_root_index);
            nodes.push(StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: current_root_index,
                    secondary_child_index: star_node_index,
                    semi_major_axis_au: outer_orbit.semi_major_axis_au,
                    eccentricity: outer_orbit.eccentricity,
                }),
            });
            current_root_index = new_root_index;
        }

        StellarHierarchy {
            root_index: current_root_index,
            nodes,
        }
    }

    fn build_hierarchy_planetary_regions(
        log_level: u8,
        stars: &mut [Star],
        hierarchy: &StellarHierarchy,
        generation_mode: GenerationMode,
    ) -> Vec<PlanetaryRegion> {
        let mut planetary_regions = Vec::new();
        for (node_index, node) in hierarchy.nodes.iter().enumerate() {
            match &node.kind {
                StellarHierarchyNodeKind::Star { star_index } => {
                    let star_ref = Self::star_ref_from_global_star_index(*star_index);
                    let star = stars[*star_index].clone();
                    let hierarchy_host = Some(Self::hierarchy_star_host(node_index, *star_index));
                    let Some(parent_context) = Self::hierarchy_parent_binary_context(stars, hierarchy, node_index)
                    else {
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR suppressing hierarchy circumstellar host node={} star={} because no parent binary context was available; leaving an empty disk to avoid inventing unsupported bounds",
                            node_index,
                            star_index
                        );
                        stars[*star_index].body.accretion_disk =
                            Some(Arc::new(RwLock::new(Self::empty_unbounded_circumstellar_disk(&star))));
                        continue;
                    };

                    log!(
                        log_level,
                        1,
                        "DEBUG-MULTISTAR discovered hierarchy circumstellar host node={} star={} hosted_system={} local_star={} parent_node={} sibling_node={}",
                        node_index,
                        star_index,
                        star_ref.hosted_system_index,
                        star_ref.star_index,
                        parent_context.parent_node_index,
                        parent_context.sibling_node_index
                    );

                    let region_result = if star_ref.hosted_system_index == 0 {
                        Self::build_circumstellar_region(
                            &star,
                            &parent_context.topology,
                            star_ref.hosted_system_index,
                            star_ref.star_index,
                            generation_mode,
                            hierarchy_host,
                        )
                    } else {
                        Self::build_outer_companion_circumstellar_region(
                            &star,
                            &parent_context.topology,
                            star_ref.hosted_system_index,
                            generation_mode,
                            hierarchy_host,
                        )
                    };

                    if let Some((region, disk)) = region_result {
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR accepted hierarchy circumstellar region node={} star={} inner={:.3} AU outer={:.3} AU stability_outer={:.3} AU truncation_outer={:.3} AU",
                            node_index,
                            star_index,
                            region.stable_region.inner_au,
                            region.stable_region.outer_au,
                            region.stable_region.stability_outer_au,
                            region.stable_region.truncation_outer_au
                        );
                        stars[*star_index].body.accretion_disk = Some(Arc::new(RwLock::new(disk)));
                        planetary_regions.push(region);
                    } else {
                        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
                        let bounds = Self::circumstellar_outer_bounds(
                            default_disk.planet_outer_bound,
                            parent_context.topology.semi_major_axis_au,
                            parent_context.topology.eccentricity,
                            parent_context.topology.mass_ratio,
                        );
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR rejected hierarchy circumstellar region node={} star={} stability_outer={:.3} AU truncation_outer={:.3} AU final_outer={:.3} AU width_below_minimum={:.3} AU",
                            node_index,
                            star_index,
                            bounds.stability_outer_au,
                            bounds.truncation_outer_au,
                            bounds.outer_au,
                            MIN_PLANETARY_REGION_WIDTH_AU
                        );
                        let disk = if star_ref.hosted_system_index == 0 {
                            Self::empty_circumstellar_disk(&star, &parent_context.topology, star_ref.star_index)
                        } else {
                            Self::empty_outer_companion_circumstellar_disk(&star, &parent_context.topology)
                        };
                        stars[*star_index].body.accretion_disk = Some(Arc::new(RwLock::new(disk)));
                    }
                }
                StellarHierarchyNodeKind::Binary(binary) => {
                    let hierarchy_host = Some(Self::hierarchy_binary_host(node_index, binary));
                    let Some(PlanetHost::Circumbinary {
                        hosted_system_index,
                        primary_star_index,
                        secondary_star_index,
                    }) = Self::hierarchy_binary_compatibility_host(hierarchy, node_index)
                    else {
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR suppressing hierarchy binary node {} circumbinary region because no compatibility PlanetHost projection exists for its canonical hierarchy host; cross-host circummultiple runtime support remains deferred",
                            node_index
                        );
                        continue;
                    };

                    let Some(inner_binary) = Self::hierarchy_binary_topology(stars, hierarchy, node_index) else {
                        continue;
                    };
                    let StellarHierarchyNodeKind::Star {
                        star_index: primary_global_star_index,
                    } = hierarchy.nodes[binary.primary_child_index].kind
                    else {
                        continue;
                    };
                    let StellarHierarchyNodeKind::Star {
                        star_index: secondary_global_star_index,
                    } = hierarchy.nodes[binary.secondary_child_index].kind
                    else {
                        continue;
                    };
                    let primary = stars[primary_global_star_index].clone();
                    let secondary = stars[secondary_global_star_index].clone();

                    log!(
                        log_level,
                        1,
                        "DEBUG-MULTISTAR discovered hierarchy circumbinary host node={} hosted_system={} stars=({}, {})",
                        node_index,
                        hosted_system_index,
                        primary_star_index,
                        secondary_star_index
                    );

                    let region = if let Some(parent_context) =
                        Self::hierarchy_parent_binary_context(stars, hierarchy, node_index)
                    {
                        Self::build_outer_limited_circumbinary_region(
                            &inner_binary,
                            &parent_context.topology,
                            hosted_system_index,
                            &primary,
                            &secondary,
                            generation_mode,
                            hierarchy_host,
                        )
                    } else {
                        Self::build_circumbinary_region(
                            &inner_binary,
                            hosted_system_index,
                            &primary,
                            &secondary,
                            generation_mode,
                            hierarchy_host,
                        )
                    };

                    if let Some(region) = region {
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR accepted hierarchy circumbinary region node={} inner={:.3} AU outer={:.3} AU stability_inner={:.3} AU truncation_inner={:.3} AU stability_outer={:.3} AU truncation_outer={:.3} AU",
                            node_index,
                            region.stable_region.inner_au,
                            region.stable_region.outer_au,
                            region.stable_region.stability_inner_au,
                            region.stable_region.truncation_inner_au,
                            region.stable_region.stability_outer_au,
                            region.stable_region.truncation_outer_au
                        );
                        planetary_regions.push(region);
                    } else {
                        let aggregate = Self::hierarchy_node_aggregate(stars, hierarchy, node_index);
                        let default_outer_au = 50.0 * aggregate.total_mass_solar.cbrt();
                        let outer_bounds = Self::hierarchy_parent_binary_context(stars, hierarchy, node_index)
                            .map(|parent_context| {
                                Self::circumstellar_outer_bounds(
                                    default_outer_au,
                                    parent_context.topology.semi_major_axis_au,
                                    parent_context.topology.eccentricity,
                                    parent_context.topology.mass_ratio,
                                )
                            })
                            .unwrap_or(CircumstellarOuterBounds {
                                outer_au: default_outer_au,
                                stability_outer_au: default_outer_au,
                                truncation_outer_au: default_outer_au,
                            });
                        log!(
                            log_level,
                            1,
                            "DEBUG-MULTISTAR rejected hierarchy circumbinary region node={} stability_inner={:.3} AU truncation_inner={:.3} AU stability_outer={:.3} AU truncation_outer={:.3} AU width_below_minimum={:.3} AU",
                            node_index,
                            Self::circumbinary_stability_inner_au(
                                inner_binary.semi_major_axis_au,
                                inner_binary.eccentricity,
                                inner_binary.mass_ratio,
                            ),
                            Self::circumbinary_truncation_inner_au(inner_binary.semi_major_axis_au, inner_binary.eccentricity),
                            outer_bounds.stability_outer_au,
                            outer_bounds.truncation_outer_au,
                            MIN_PLANETARY_REGION_WIDTH_AU
                        );
                    }
                }
                StellarHierarchyNodeKind::Unresolved { .. } => {
                    log!(
                        log_level,
                        1,
                        "DEBUG-MULTISTAR suppressing unresolved hierarchy node {} because Slice 3 only emits regions for resolved star and direct-binary hosts",
                        node_index
                    );
                }
            }
        }

        Self::sort_planetary_regions(&mut planetary_regions);
        planetary_regions
    }

    fn hierarchy_node_aggregate(stars: &[Star], hierarchy: &StellarHierarchy, node_index: usize) -> HierarchyAggregate {
        match &hierarchy.nodes[node_index].kind {
            StellarHierarchyNodeKind::Star { star_index } => {
                let star = &stars[*star_index];
                HierarchyAggregate {
                    total_mass_solar: star.mass_in_sols,
                    total_luminosity_solar: star.luminosity_in_sols,
                }
            }
            StellarHierarchyNodeKind::Binary(binary) => {
                let primary = Self::hierarchy_node_aggregate(stars, hierarchy, binary.primary_child_index);
                let secondary = Self::hierarchy_node_aggregate(stars, hierarchy, binary.secondary_child_index);
                HierarchyAggregate {
                    total_mass_solar: primary.total_mass_solar + secondary.total_mass_solar,
                    total_luminosity_solar: primary.total_luminosity_solar + secondary.total_luminosity_solar,
                }
            }
            StellarHierarchyNodeKind::Unresolved { star_indices } => star_indices.iter().fold(
                HierarchyAggregate {
                    total_mass_solar: 0.0,
                    total_luminosity_solar: 0.0,
                },
                |aggregate, star_index| HierarchyAggregate {
                    total_mass_solar: aggregate.total_mass_solar + stars[*star_index].mass_in_sols,
                    total_luminosity_solar: aggregate.total_luminosity_solar + stars[*star_index].luminosity_in_sols,
                },
            ),
        }
    }

    fn hierarchy_parent_binary_context(
        stars: &[Star],
        hierarchy: &StellarHierarchy,
        node_index: usize,
    ) -> Option<HierarchyParentBinaryContext> {
        let parent_node_index = hierarchy.nodes.get(node_index)?.parent_index?;
        let StellarHierarchyNodeKind::Binary(parent_binary) = &hierarchy.nodes.get(parent_node_index)?.kind else {
            return None;
        };
        let sibling_node_index = if parent_binary.primary_child_index == node_index {
            parent_binary.secondary_child_index
        } else if parent_binary.secondary_child_index == node_index {
            parent_binary.primary_child_index
        } else {
            return None;
        };

        let host_aggregate = Self::hierarchy_node_aggregate(stars, hierarchy, node_index);
        let sibling_aggregate = Self::hierarchy_node_aggregate(stars, hierarchy, sibling_node_index);
        let total_mass_solar = host_aggregate.total_mass_solar + sibling_aggregate.total_mass_solar;
        Some(HierarchyParentBinaryContext {
            parent_node_index,
            sibling_node_index,
            topology: BinaryTopology {
                primary_index: 0,
                secondary_index: 1,
                semi_major_axis_au: parent_binary.semi_major_axis_au,
                eccentricity: parent_binary.eccentricity,
                total_mass_solar,
                mass_ratio: if total_mass_solar > 0.0 {
                    sibling_aggregate.total_mass_solar / total_mass_solar
                } else {
                    0.0
                },
            },
        })
    }

    fn hierarchy_binary_topology(
        stars: &[Star],
        hierarchy: &StellarHierarchy,
        node_index: usize,
    ) -> Option<BinaryTopology> {
        let StellarHierarchyNodeKind::Binary(binary) = &hierarchy.nodes.get(node_index)?.kind else {
            return None;
        };
        let primary_aggregate = Self::hierarchy_node_aggregate(stars, hierarchy, binary.primary_child_index);
        let secondary_aggregate = Self::hierarchy_node_aggregate(stars, hierarchy, binary.secondary_child_index);
        let total_mass_solar = primary_aggregate.total_mass_solar + secondary_aggregate.total_mass_solar;

        Some(BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: binary.semi_major_axis_au,
            eccentricity: binary.eccentricity,
            total_mass_solar,
            mass_ratio: if total_mass_solar > 0.0 {
                secondary_aggregate.total_mass_solar / total_mass_solar
            } else {
                0.0
            },
        })
    }

    fn hierarchy_binary_compatibility_host(hierarchy: &StellarHierarchy, node_index: usize) -> Option<PlanetHost> {
        let StellarHierarchyNodeKind::Binary(binary) = &hierarchy.nodes.get(node_index)?.kind else {
            return None;
        };
        let StellarHierarchyNodeKind::Star {
            star_index: primary_global_star_index,
        } = hierarchy.nodes.get(binary.primary_child_index)?.kind
        else {
            return None;
        };
        let StellarHierarchyNodeKind::Star {
            star_index: secondary_global_star_index,
        } = hierarchy.nodes.get(binary.secondary_child_index)?.kind
        else {
            return None;
        };

        let primary_star_ref = Self::star_ref_from_global_star_index(primary_global_star_index);
        let secondary_star_ref = Self::star_ref_from_global_star_index(secondary_global_star_index);
        if primary_star_ref.hosted_system_index != secondary_star_ref.hosted_system_index {
            return None;
        }

        Some(PlanetHost::Circumbinary {
            hosted_system_index: primary_star_ref.hosted_system_index,
            primary_star_index: primary_star_ref.star_index,
            secondary_star_index: secondary_star_ref.star_index,
        })
    }

    fn star_ref_from_global_star_index(global_star_index: usize) -> StarRef {
        if global_star_index < 2 {
            StarRef {
                hosted_system_index: 0,
                star_index: global_star_index,
            }
        } else {
            StarRef {
                hosted_system_index: global_star_index - 1,
                star_index: 0,
            }
        }
    }

    fn sort_planetary_regions(planetary_regions: &mut [PlanetaryRegion]) {
        planetary_regions.sort_by(|left, right| {
            let left_key = left.stable_region.host.display_order_key();
            let right_key = right.stable_region.host.display_order_key();
            left_key.cmp(&right_key).then_with(|| {
                left.stable_region
                    .inner_au
                    .partial_cmp(&right.stable_region.inner_au)
                    .unwrap()
            })
        });
    }

    fn build_binary_planetary_regions(
        stars: &mut [Star],
        topology: &StellarTopology,
        hosted_system_index: usize,
        generation_mode: GenerationMode,
    ) -> Vec<PlanetaryRegion> {
        let StellarTopology::Binary(binary) = topology else {
            return Vec::new();
        };

        let mut planetary_regions = Vec::new();
        let primary = stars[0].clone();
        let secondary = stars[1].clone();

        if let Some((region, disk)) = Self::build_circumstellar_region(
            &primary,
            binary,
            hosted_system_index,
            0,
            generation_mode,
            Some(HierarchyHostRef::Star {
                node_index: 0,
                global_star_index: 0,
            }),
        ) {
            stars[0].body.accretion_disk = Some(Arc::new(RwLock::new(disk)));
            planetary_regions.push(region);
        } else {
            stars[0].body.accretion_disk = Some(Arc::new(RwLock::new(Self::empty_circumstellar_disk(
                &primary, binary, 0,
            ))));
        }

        if let Some((region, disk)) = Self::build_circumstellar_region(
            &secondary,
            binary,
            hosted_system_index,
            1,
            generation_mode,
            Some(HierarchyHostRef::Star {
                node_index: 1,
                global_star_index: 1,
            }),
        ) {
            stars[1].body.accretion_disk = Some(Arc::new(RwLock::new(disk)));
            planetary_regions.push(region);
        } else {
            stars[1].body.accretion_disk = Some(Arc::new(RwLock::new(Self::empty_circumstellar_disk(
                &secondary, binary, 1,
            ))));
        }

        if let Some(region) = Self::build_circumbinary_region(
            binary,
            hosted_system_index,
            &primary,
            &secondary,
            generation_mode,
            Some(HierarchyHostRef::Binary {
                node_index: 2,
                primary_child_index: 0,
                secondary_child_index: 1,
            }),
        ) {
            planetary_regions.push(region);
        }

        Self::sort_planetary_regions(&mut planetary_regions);
        planetary_regions
    }

    fn configure_star_orbit(star: &mut Star, a: f64, e: f64) {
        star.a = a;
        star.body.a = a;
        star.body.e = e;
    }

    fn build_circumstellar_region(
        star: &Star,
        binary: &BinaryTopology,
        hosted_system_index: usize,
        star_index: usize,
        generation_mode: GenerationMode,
        hierarchy_host: Option<HierarchyHostRef>,
    ) -> Option<(PlanetaryRegion, AccretionDisk)> {
        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
        let companion_mass = binary.total_mass_solar - star.mass_in_sols;
        let companion_fraction = if binary.total_mass_solar > 0.0 {
            companion_mass / binary.total_mass_solar
        } else {
            0.0
        };
        let inner_au = default_disk.planet_inner_bound;
        let outer_bounds = Self::circumstellar_outer_bounds(
            default_disk.planet_outer_bound,
            binary.semi_major_axis_au,
            binary.eccentricity,
            companion_fraction,
        );
        let outer_au = outer_bounds.outer_au;
        let disk = AccretionDisk::new_with_planet_bounds(
            star.mass_in_sols,
            star.luminosity_in_sols,
            MassType::Star,
            inner_au,
            outer_au,
        );

        if outer_au - inner_au < MIN_PLANETARY_REGION_WIDTH_AU {
            return None;
        }

        let mut bounded_disk = disk;
        bounded_disk.accrete_with_mode(generation_mode);
        let region = Self::build_planetary_region(
            StableRegion {
                host: PlanetHost::Circumstellar {
                    star_ref: StarRef {
                        hosted_system_index,
                        star_index,
                    },
                },
                hierarchy_host,
                inner_au,
                outer_au,
                stability_inner_au: default_disk.planet_inner_bound,
                stability_outer_au: outer_bounds.stability_outer_au,
                truncation_inner_au: default_disk.planet_inner_bound,
                truncation_outer_au: outer_bounds.truncation_outer_au,
            },
            star.mass_in_sols,
            star.luminosity_in_sols,
            bounded_disk.bodies.clone(),
        );

        Some((region, bounded_disk))
    }

    fn empty_circumstellar_disk(star: &Star, binary: &BinaryTopology, star_index: usize) -> AccretionDisk {
        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
        let companion_mass = binary.total_mass_solar - star.mass_in_sols;
        let companion_fraction = if binary.total_mass_solar > 0.0 {
            companion_mass / binary.total_mass_solar
        } else {
            0.0
        };
        let outer_au = Self::circumstellar_outer_bounds(
            default_disk.planet_outer_bound,
            binary.semi_major_axis_au,
            binary.eccentricity,
            companion_fraction,
        )
        .outer_au;
        let _ = star_index;
        AccretionDisk::new_with_planet_bounds(
            star.mass_in_sols,
            star.luminosity_in_sols,
            MassType::Star,
            default_disk.planet_inner_bound,
            outer_au,
        )
    }

    fn build_outer_companion_circumstellar_region(
        star: &Star,
        outer_binary: &BinaryTopology,
        hosted_system_index: usize,
        generation_mode: GenerationMode,
        hierarchy_host: Option<HierarchyHostRef>,
    ) -> Option<(PlanetaryRegion, AccretionDisk)> {
        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
        let companion_mass = outer_binary.total_mass_solar - star.mass_in_sols;
        let companion_fraction = if outer_binary.total_mass_solar > 0.0 {
            companion_mass / outer_binary.total_mass_solar
        } else {
            0.0
        };
        let outer_bounds = Self::circumstellar_outer_bounds(
            default_disk.planet_outer_bound,
            outer_binary.semi_major_axis_au,
            outer_binary.eccentricity,
            companion_fraction,
        );
        let inner_au = default_disk.planet_inner_bound;
        let outer_au = outer_bounds.outer_au;
        if outer_au - inner_au < MIN_PLANETARY_REGION_WIDTH_AU {
            return None;
        }

        let mut bounded_disk = AccretionDisk::new_with_planet_bounds(
            star.mass_in_sols,
            star.luminosity_in_sols,
            MassType::Star,
            inner_au,
            outer_au,
        );
        bounded_disk.accrete_with_mode(generation_mode);
        let region = Self::build_planetary_region(
            StableRegion {
                host: PlanetHost::Circumstellar {
                    star_ref: StarRef {
                        hosted_system_index,
                        star_index: 0,
                    },
                },
                hierarchy_host,
                inner_au,
                outer_au,
                stability_inner_au: default_disk.planet_inner_bound,
                stability_outer_au: outer_bounds.stability_outer_au,
                truncation_inner_au: default_disk.planet_inner_bound,
                truncation_outer_au: outer_bounds.truncation_outer_au,
            },
            star.mass_in_sols,
            star.luminosity_in_sols,
            bounded_disk.bodies.clone(),
        );

        Some((region, bounded_disk))
    }

    fn empty_outer_companion_circumstellar_disk(star: &Star, outer_binary: &BinaryTopology) -> AccretionDisk {
        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
        let companion_mass = outer_binary.total_mass_solar - star.mass_in_sols;
        let companion_fraction = if outer_binary.total_mass_solar > 0.0 {
            companion_mass / outer_binary.total_mass_solar
        } else {
            0.0
        };
        let outer_au = Self::circumstellar_outer_bounds(
            default_disk.planet_outer_bound,
            outer_binary.semi_major_axis_au,
            outer_binary.eccentricity,
            companion_fraction,
        )
        .outer_au;

        AccretionDisk::new_with_planet_bounds(
            star.mass_in_sols,
            star.luminosity_in_sols,
            MassType::Star,
            default_disk.planet_inner_bound,
            outer_au,
        )
    }

    fn empty_unbounded_circumstellar_disk(star: &Star) -> AccretionDisk {
        let default_disk = AccretionDisk::new(&star.body, star.luminosity_in_sols);
        AccretionDisk::new_with_planet_bounds(
            star.mass_in_sols,
            star.luminosity_in_sols,
            MassType::Star,
            default_disk.planet_inner_bound,
            default_disk.planet_inner_bound,
        )
    }

    fn build_circumbinary_region(
        binary: &BinaryTopology,
        hosted_system_index: usize,
        primary: &Star,
        secondary: &Star,
        generation_mode: GenerationMode,
        hierarchy_host: Option<HierarchyHostRef>,
    ) -> Option<PlanetaryRegion> {
        let host_mass_solar = primary.mass_in_sols + secondary.mass_in_sols;
        let host_luminosity_solar = primary.luminosity_in_sols + secondary.luminosity_in_sols;
        let default_inner_au = 0.3 * host_mass_solar.cbrt();
        let default_outer_au = 50.0 * host_mass_solar.cbrt();
        let stability_inner_au =
            Self::circumbinary_stability_inner_au(binary.semi_major_axis_au, binary.eccentricity, binary.mass_ratio);
        let truncation_inner_au =
            Self::circumbinary_truncation_inner_au(binary.semi_major_axis_au, binary.eccentricity);
        let inner_au = default_inner_au.max(stability_inner_au).max(truncation_inner_au);
        let outer_au = default_outer_au;
        if outer_au - inner_au < MIN_PLANETARY_REGION_WIDTH_AU {
            return None;
        }

        let mut bounded_disk = AccretionDisk::new_with_planet_bounds(
            host_mass_solar,
            host_luminosity_solar,
            MassType::Star,
            inner_au,
            outer_au,
        );
        bounded_disk.accrete_with_mode(generation_mode);

        Some(Self::build_planetary_region(
            StableRegion {
                host: PlanetHost::Circumbinary {
                    hosted_system_index,
                    primary_star_index: binary.primary_index,
                    secondary_star_index: binary.secondary_index,
                },
                hierarchy_host,
                inner_au,
                outer_au,
                stability_inner_au,
                stability_outer_au: default_outer_au,
                truncation_inner_au,
                truncation_outer_au: default_outer_au,
            },
            host_mass_solar,
            host_luminosity_solar,
            bounded_disk.bodies.clone(),
        ))
    }

    fn build_outer_limited_circumbinary_region(
        inner_binary: &BinaryTopology,
        outer_binary: &BinaryTopology,
        hosted_system_index: usize,
        primary: &Star,
        secondary: &Star,
        generation_mode: GenerationMode,
        hierarchy_host: Option<HierarchyHostRef>,
    ) -> Option<PlanetaryRegion> {
        let host_mass_solar = primary.mass_in_sols + secondary.mass_in_sols;
        let host_luminosity_solar = primary.luminosity_in_sols + secondary.luminosity_in_sols;
        let default_inner_au = 0.3 * host_mass_solar.cbrt();
        let default_outer_au = 50.0 * host_mass_solar.cbrt();
        let stability_inner_au = Self::circumbinary_stability_inner_au(
            inner_binary.semi_major_axis_au,
            inner_binary.eccentricity,
            inner_binary.mass_ratio,
        );
        let truncation_inner_au =
            Self::circumbinary_truncation_inner_au(inner_binary.semi_major_axis_au, inner_binary.eccentricity);
        let outer_bounds = Self::circumstellar_outer_bounds(
            default_outer_au,
            outer_binary.semi_major_axis_au,
            outer_binary.eccentricity,
            outer_binary.mass_ratio,
        );
        let inner_au = default_inner_au.max(stability_inner_au).max(truncation_inner_au);
        let outer_au = outer_bounds.outer_au;
        if outer_au - inner_au < MIN_PLANETARY_REGION_WIDTH_AU {
            return None;
        }

        let mut bounded_disk = AccretionDisk::new_with_planet_bounds(
            host_mass_solar,
            host_luminosity_solar,
            MassType::Star,
            inner_au,
            outer_au,
        );
        bounded_disk.accrete_with_mode(generation_mode);

        Some(Self::build_planetary_region(
            StableRegion {
                host: PlanetHost::Circumbinary {
                    hosted_system_index,
                    primary_star_index: inner_binary.primary_index,
                    secondary_star_index: inner_binary.secondary_index,
                },
                hierarchy_host,
                inner_au,
                outer_au,
                stability_inner_au,
                stability_outer_au: outer_bounds.stability_outer_au,
                truncation_inner_au,
                truncation_outer_au: outer_bounds.truncation_outer_au,
            },
            host_mass_solar,
            host_luminosity_solar,
            bounded_disk.bodies.clone(),
        ))
    }

    fn build_planetary_region(
        stable_region: StableRegion,
        host_mass_solar: f64,
        host_luminosity_solar: f64,
        mut bodies: Vec<crate::body::Body>,
    ) -> PlanetaryRegion {
        bodies.sort_by(|left, right| left.a.partial_cmp(&right.a).unwrap_or(std::cmp::Ordering::Equal));
        let bodies = Self::apply_spacing_cleanup(host_mass_solar, bodies);
        PlanetaryRegion::new(stable_region, host_mass_solar, host_luminosity_solar, bodies)
    }

    fn apply_spacing_cleanup(host_mass_solar: f64, bodies: Vec<crate::body::Body>) -> Vec<crate::body::Body> {
        let mut kept: Vec<(usize, crate::body::Body)> = Vec::with_capacity(bodies.len());

        for (sorted_index, body) in bodies.into_iter().enumerate() {
            let candidate_index = sorted_index;
            let candidate = body;

            loop {
                let Some((last_index, last_body)) = kept.last() else {
                    kept.push((candidate_index, candidate));
                    break;
                };

                if !Self::bodies_are_overpacked(host_mass_solar, last_body, &candidate) {
                    kept.push((candidate_index, candidate));
                    break;
                }

                if Self::spacing_cleanup_prefers_outer(*last_index, last_body, candidate_index, &candidate) {
                    kept.pop();
                    continue;
                }

                break;
            }
        }

        kept.into_iter().map(|(_, body)| body).collect()
    }

    fn bodies_are_overpacked(host_mass_solar: f64, inner: &crate::body::Body, outer: &crate::body::Body) -> bool {
        let ratio = if inner.a > 0.0 {
            outer.a / inner.a
        } else {
            f64::INFINITY
        };
        if ratio < MIN_PLANETARY_LOG_SPACING_RATIO {
            return true;
        }

        let combined_mass_solar = (inner.mass_in_sols + outer.mass_in_sols).max(0.0);
        if host_mass_solar <= 0.0 || combined_mass_solar <= 0.0 {
            return false;
        }

        let average_a = ((inner.a + outer.a) * 0.5).max(0.0);
        let mutual_hill_radius = average_a * (combined_mass_solar / (3.0 * host_mass_solar)).cbrt();
        let required_spacing = MIN_PLANETARY_MUTUAL_HILL_SPACING * mutual_hill_radius;
        let actual_spacing = (outer.a - inner.a).max(0.0);

        actual_spacing < required_spacing
    }

    fn spacing_cleanup_prefers_outer(
        inner_index: usize,
        inner: &crate::body::Body,
        outer_index: usize,
        outer: &crate::body::Body,
    ) -> bool {
        let mass_scale = inner.mass_in_sols.abs().max(outer.mass_in_sols.abs()).max(1.0);
        let mass_difference = outer.mass_in_sols - inner.mass_in_sols;
        if mass_difference.abs() > PLANETARY_MASS_TIE_RELATIVE_TOLERANCE * mass_scale {
            return mass_difference > 0.0;
        }

        if outer.a > inner.a {
            return false;
        }

        if outer.a < inner.a {
            return true;
        }

        outer_index < inner_index
    }

    fn circumstellar_outer_bounds(
        default_outer_au: f64,
        binary_separation_au: f64,
        binary_eccentricity: f64,
        companion_fraction: f64,
    ) -> CircumstellarOuterBounds {
        let stability_outer_au =
            Self::circumstellar_stability_outer_au(binary_separation_au, binary_eccentricity, companion_fraction);
        let truncation_outer_au = Self::circumstellar_truncation_outer_au(binary_separation_au, binary_eccentricity);

        CircumstellarOuterBounds {
            outer_au: default_outer_au.min(stability_outer_au).min(truncation_outer_au),
            stability_outer_au,
            truncation_outer_au,
        }
    }

    fn circumstellar_stability_outer_au(
        binary_separation_au: f64,
        binary_eccentricity: f64,
        companion_fraction: f64,
    ) -> f64 {
        let e = binary_eccentricity.clamp(0.0, 0.8);
        let mu = companion_fraction.clamp(0.0, 0.95);
        let ratio = 0.464 - 0.38 * mu - 0.631 * e + 0.586 * mu * e + 0.15 * e * e - 0.198 * mu * e * e;
        (binary_separation_au * ratio).max(0.0)
    }

    fn circumbinary_stability_inner_au(binary_separation_au: f64, binary_eccentricity: f64, mass_ratio: f64) -> f64 {
        let e = binary_eccentricity.clamp(0.0, 0.8);
        let mu = mass_ratio.clamp(0.0, 0.5);
        let ratio = 1.6 + 5.1 * e - 2.22 * e * e + 4.12 * mu - 4.27 * e * mu - 5.09 * mu * mu + 4.61 * e * e * mu * mu;
        (binary_separation_au * ratio).max(0.0)
    }

    fn circumstellar_truncation_outer_au(binary_separation_au: f64, binary_eccentricity: f64) -> f64 {
        (0.3 * binary_separation_au * (1.0 - binary_eccentricity).max(0.0)).max(0.0)
    }

    fn circumbinary_truncation_inner_au(binary_separation_au: f64, binary_eccentricity: f64) -> f64 {
        (2.0 * binary_separation_au * (1.0 + binary_eccentricity.max(0.0))).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Body;
    use crate::condensation::CondensationThresholds;
    use crate::star::{LuminosityClass, SpectralClass};

    fn hierarchy(system: &StarSystem) -> &StellarHierarchy {
        system.hierarchy.as_ref().expect("expected canonical hierarchy")
    }

    fn test_body(semi_major_axis_au: f64, mass_in_sols: f64, mass_type: MassType) -> Body {
        Body::new(
            semi_major_axis_au,
            0.0,
            mass_in_sols,
            mass_type,
            6_400.0,
            1.0,
            1.0,
            None,
        )
    }

    fn test_stable_region(host: PlanetHost) -> StableRegion {
        StableRegion {
            host,
            hierarchy_host: None,
            inner_au: 0.1,
            outer_au: 20.0,
            stability_inner_au: 0.1,
            stability_outer_au: 20.0,
            truncation_inner_au: 0.1,
            truncation_outer_au: 20.0,
        }
    }

    #[test]
    fn from_stars_builds_single_topology_and_legacy_region() {
        let mut star = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        star.accrete();

        let system = StarSystem::from_stars(vec![star]);
        assert!(matches!(system.topology, StellarTopology::Single { primary_index: 0 }));
        assert_eq!(system.planetary_regions.len(), 1);
        assert!(matches!(
            system.planetary_regions[0].stable_region.host,
            PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 0,
                },
            }
        ));
        assert_eq!(
            system.planetary_regions[0].stable_region.hierarchy_host,
            Some(HierarchyHostRef::Star {
                node_index: 0,
                global_star_index: 0,
            })
        );
        assert_eq!(system.hierarchy, Some(StellarHierarchy::single(0)));
    }

    #[test]
    fn circumstellar_region_thermal_profile_uses_host_luminosity() {
        let mut star = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        star.accrete();

        let system = StarSystem::from_stars(vec![star]);
        let region = &system.planetary_regions[0];

        assert_eq!(
            region.thermal_profile,
            PlanetaryRegion::derive_thermal_profile(
                &region.stable_region,
                region.host_mass_solar,
                region.host_luminosity_solar,
            )
        );
        assert_eq!(
            region.thermal_profile.host_luminosity_solar,
            region.host_luminosity_solar
        );
    }

    #[test]
    fn build_planetary_region_preserves_well_spaced_bodies() {
        let region = StarSystem::build_planetary_region(
            test_stable_region(PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 0,
                },
            }),
            1.0,
            1.0,
            vec![
                test_body(0.5, 1.0e-6, MassType::Planet),
                test_body(0.9, 1.2e-6, MassType::Planet),
                test_body(1.5, 1.1e-6, MassType::Planet),
            ],
        );

        assert_eq!(
            region.bodies.iter().map(|body| body.a).collect::<Vec<f64>>(),
            vec![0.5, 0.9, 1.5]
        );
    }

    #[test]
    fn build_planetary_region_keeps_inner_body_when_close_pair_is_mass_tied() {
        let region = StarSystem::build_planetary_region(
            test_stable_region(PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 0,
                },
            }),
            1.0,
            1.0,
            vec![
                test_body(1.0, 1.0e-6, MassType::Planet),
                test_body(1.08, 1.0e-6, MassType::Planet),
                test_body(1.6, 1.0e-6, MassType::Planet),
            ],
        );

        assert_eq!(
            region.bodies.iter().map(|body| body.a).collect::<Vec<f64>>(),
            vec![1.0, 1.6]
        );
    }

    #[test]
    fn build_planetary_region_keeps_more_massive_outer_body_when_close_pair_is_overpacked() {
        let region = StarSystem::build_planetary_region(
            test_stable_region(PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 0,
                },
            }),
            1.0,
            1.0,
            vec![
                test_body(1.0, 1.0e-6, MassType::Planet),
                test_body(1.08, 3.0e-6, MassType::Planet),
                test_body(1.8, 1.0e-6, MassType::Planet),
            ],
        );

        assert_eq!(
            region.bodies.iter().map(|body| body.a).collect::<Vec<f64>>(),
            vec![1.08, 1.8]
        );
    }

    #[test]
    fn build_planetary_region_uses_mutual_hill_proxy_for_massive_pairs() {
        let region = StarSystem::build_planetary_region(
            test_stable_region(PlanetHost::Circumbinary {
                hosted_system_index: 0,
                primary_star_index: 0,
                secondary_star_index: 1,
            }),
            1.0,
            1.0,
            vec![
                test_body(5.0, 9.0e-4, MassType::GasGiant),
                test_body(6.0, 1.1e-3, MassType::GasGiant),
                test_body(10.0, 8.0e-4, MassType::GasGiant),
            ],
        );

        assert_eq!(
            region.bodies.iter().map(|body| body.a).collect::<Vec<f64>>(),
            vec![6.0, 10.0]
        );
        assert!(matches!(region.stable_region.host, PlanetHost::Circumbinary { .. }));
    }

    #[test]
    fn binary_topology_uses_descending_mass_order() {
        let mut primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 0, 0.0);
        let mut secondary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 5, 0.0);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 12.0, 0.2);

        let topology = StarSystem::binary_topology_from_stars(&[primary.clone(), secondary.clone()]);
        assert_eq!(topology.primary_index, 0);
        assert_eq!(topology.secondary_index, 1);
        assert!(primary.mass_in_sols >= secondary.mass_in_sols);
        assert_eq!(topology.semi_major_axis_au, 12.0);
    }

    #[test]
    fn from_stars_builds_binary_hierarchy_metadata() {
        let mut primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let mut secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 12.0);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 12.0, 0.2);

        let system = StarSystem::from_stars(vec![primary, secondary]);
        let hierarchy = hierarchy(&system);

        assert_eq!(hierarchy.root_index, 2);
        assert_eq!(hierarchy.nodes.len(), 3);
        assert!(matches!(
            hierarchy.nodes[0],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 0 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[1],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 1 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[2],
            StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 0,
                    secondary_child_index: 1,
                    semi_major_axis_au: 12.0,
                    eccentricity: 0.2,
                }),
            }
        ));
        assert!(system.hosted_systems[0]
            .planetary_regions
            .iter()
            .all(|region| region.stable_region.hierarchy_host.is_some()));
    }

    #[test]
    fn circumbinary_region_thermal_profile_uses_combined_host_luminosity() {
        let primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 0.0);
        let binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: 8.0,
            eccentricity: 0.2,
            total_mass_solar: primary.mass_in_sols + secondary.mass_in_sols,
            mass_ratio: secondary.mass_in_sols / (primary.mass_in_sols + secondary.mass_in_sols),
        };

        let region =
            StarSystem::build_circumbinary_region(&binary, 0, &primary, &secondary, GenerationMode::Aggregation, None)
                .expect("expected moderate binary to produce a circumbinary region");
        assert!(region.stable_region.inner_au < region.stable_region.outer_au);
        assert!(matches!(region.stable_region.host, PlanetHost::Circumbinary { .. }));
        assert_eq!(
            region.thermal_profile,
            PlanetaryRegion::derive_thermal_profile(
                &region.stable_region,
                region.host_mass_solar,
                region.host_luminosity_solar,
            )
        );
        assert_eq!(
            region.host_luminosity_solar,
            primary.luminosity_in_sols + secondary.luminosity_in_sols
        );
        assert_eq!(
            region.thermal_profile.host_luminosity_solar,
            primary.luminosity_in_sols + secondary.luminosity_in_sols
        );
    }

    #[test]
    fn circumbinary_sampled_profile_changes_with_host_luminosity() {
        let dim_primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let dim_secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 0.0);
        let mut bright_primary = dim_primary.clone();
        let mut bright_secondary = dim_secondary.clone();
        bright_primary.luminosity_in_sols *= 16.0;
        bright_secondary.luminosity_in_sols *= 16.0;

        let binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: 8.0,
            eccentricity: 0.2,
            total_mass_solar: dim_primary.mass_in_sols + dim_secondary.mass_in_sols,
            mass_ratio: dim_secondary.mass_in_sols / (dim_primary.mass_in_sols + dim_secondary.mass_in_sols),
        };
        let thresholds = CondensationThresholds::default();

        let dim_region = StarSystem::build_circumbinary_region(
            &binary,
            0,
            &dim_primary,
            &dim_secondary,
            GenerationMode::Aggregation,
            None,
        )
        .expect("expected dim binary to produce a circumbinary region");
        let bright_region = StarSystem::build_circumbinary_region(
            &binary,
            0,
            &bright_primary,
            &bright_secondary,
            GenerationMode::Aggregation,
            None,
        )
        .expect("expected bright binary to produce a circumbinary region");

        assert_eq!(dim_region.stable_region.inner_au, bright_region.stable_region.inner_au);
        assert_eq!(dim_region.stable_region.outer_au, bright_region.stable_region.outer_au);

        let sample_radius_au = dim_region.stable_region.inner_au + 5.0;
        let dim_inputs = dim_region
            .thermal_profile
            .sample_planet_formation_inputs(sample_radius_au, thresholds);
        let bright_inputs = bright_region
            .thermal_profile
            .sample_planet_formation_inputs(sample_radius_au, thresholds);

        assert!(bright_inputs.temperature_k > dim_inputs.temperature_k);
        assert!(bright_inputs.condensation_fractions.gas > dim_inputs.condensation_fractions.gas);
        assert!(bright_inputs.condensation_fractions.volatile_ices < dim_inputs.condensation_fractions.volatile_ices);
    }

    #[test]
    fn display_order_prefers_primary_then_secondary_then_circumbinary() {
        assert!(
            PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 0,
                },
            }
            .display_order_key()
                < PlanetHost::Circumstellar {
                    star_ref: StarRef {
                        hosted_system_index: 0,
                        star_index: 1,
                    },
                }
                .display_order_key()
        );
        assert!(
            PlanetHost::Circumstellar {
                star_ref: StarRef {
                    hosted_system_index: 0,
                    star_index: 1,
                },
            }
            .display_order_key()
                < PlanetHost::Circumbinary {
                    hosted_system_index: 0,
                    primary_star_index: 0,
                    secondary_star_index: 1,
                }
                .display_order_key()
        );
    }

    #[test]
    fn outer_companion_sampling_caps_cleanly_near_maximum_separation() {
        assert_eq!(StarSystem::sample_outer_companion_separation_au(3_000.0), 5_000.0);
        assert_eq!(StarSystem::sample_outer_companion_separation_au(5_000.0), 5_000.0);
    }

    #[test]
    fn hierarchical_triple_generation_creates_binary_and_tertiary_hosts() {
        let mut primary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 4, 0.0);
        let mut secondary = Star::new(SpectralClass::G, LuminosityClass::WhiteDwarf, 6, 34.52);
        let mut tertiary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 5, 425.78);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 34.52, 0.378);
        StarSystem::configure_star_orbit(&mut tertiary, 425.78, 0.001);

        let system =
            StarSystem::generate_triple_star_system(0, vec![primary, secondary, tertiary], GenerationMode::Aggregation);

        assert_eq!(system.hosted_systems.len(), 2);
        assert!(matches!(system.hosted_systems[0].topology, StellarTopology::Binary(_)));
        assert!(matches!(
            system.hosted_systems[1].topology,
            StellarTopology::Single { primary_index: 0 }
        ));
        assert_eq!(system.hosted_systems[1].stars[0].a, 425.78);
        assert_eq!(system.planetary_regions.len(), 3);

        let inner_regions = &system.hosted_systems[0].planetary_regions;
        let primary_region = inner_regions
            .iter()
            .find(|region| {
                matches!(
                    region.stable_region.host,
                    PlanetHost::Circumstellar {
                        star_ref: StarRef {
                            hosted_system_index: 0,
                            star_index: 0,
                        },
                    }
                )
            })
            .expect("expected inner primary region");
        let secondary_region = inner_regions
            .iter()
            .find(|region| {
                matches!(
                    region.stable_region.host,
                    PlanetHost::Circumstellar {
                        star_ref: StarRef {
                            hosted_system_index: 0,
                            star_index: 1,
                        },
                    }
                )
            })
            .expect("expected inner secondary region");
        let tertiary_region = &system.hosted_systems[1].planetary_regions[0];

        assert!(primary_region.stable_region.outer_au < 10.0);
        assert!(secondary_region.stable_region.outer_au < 10.0);
        assert_eq!(
            primary_region.stable_region.hierarchy_host,
            Some(HierarchyHostRef::Star {
                node_index: 0,
                global_star_index: 0,
            })
        );
        assert_eq!(
            secondary_region.stable_region.hierarchy_host,
            Some(HierarchyHostRef::Star {
                node_index: 1,
                global_star_index: 1,
            })
        );
        assert!(inner_regions.iter().all(|region| {
            !matches!(
                region.stable_region.host,
                PlanetHost::Circumbinary {
                    hosted_system_index: 0,
                    ..
                }
            )
        }));
        assert!(tertiary_region.stable_region.outer_au < 100.0);
        assert_eq!(
            tertiary_region.stable_region.hierarchy_host,
            Some(HierarchyHostRef::Star {
                node_index: 3,
                global_star_index: 2,
            })
        );

        let hierarchy = hierarchy(&system);
        assert_eq!(hierarchy.root_index, 4);
        assert_eq!(hierarchy.nodes.len(), 5);
        assert!(matches!(
            hierarchy.nodes[0],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 0 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[1],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 1 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[2],
            StellarHierarchyNode {
                parent_index: Some(4),
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 0,
                    secondary_child_index: 1,
                    semi_major_axis_au: 34.52,
                    eccentricity: 0.378,
                }),
            }
        ));
        assert!(matches!(
            hierarchy.nodes[3],
            StellarHierarchyNode {
                parent_index: Some(4),
                kind: StellarHierarchyNodeKind::Star { star_index: 2 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[4],
            StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 2,
                    secondary_child_index: 3,
                    semi_major_axis_au: 425.78,
                    eccentricity: 0.001,
                }),
            }
        ));
    }

    #[test]
    fn ambiguous_triple_generation_returns_no_planetary_regions() {
        let mut primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let mut secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 4, 30.0);
        let mut tertiary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 3, 60.0);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 30.0, 0.6);
        StarSystem::configure_star_orbit(&mut tertiary, 60.0, 0.0);

        let system =
            StarSystem::generate_triple_star_system(0, vec![primary, secondary, tertiary], GenerationMode::Aggregation);

        assert_eq!(system.hosted_systems.len(), 1);
        assert!(matches!(
            system.hosted_systems[0].topology,
            StellarTopology::HigherMultiplicity {
                primary_index: 0,
                star_count: 3,
            }
        ));
        assert!(system.planetary_regions.is_empty());
        assert!(system.hosted_systems[0].planetary_regions.is_empty());
        assert_eq!(system.hierarchy, Some(StellarHierarchy::unresolved(vec![0, 1, 2])));
    }

    #[test]
    fn resolved_four_star_generation_emits_conservative_hierarchy_justified_regions() {
        let mut primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let mut secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 10.0);
        let mut tertiary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 4, 120.0);
        let mut quaternary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 6, 1_200.0);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 10.0, 0.1);
        StarSystem::configure_star_orbit(&mut tertiary, 120.0, 0.05);
        StarSystem::configure_star_orbit(&mut quaternary, 1_200.0, 0.02);

        let system = StarSystem::generate_higher_multiplicity_star_system(
            0,
            vec![primary, secondary, tertiary, quaternary],
            GenerationMode::Aggregation,
        );

        assert_eq!(system.hosted_systems.len(), 3);
        assert!(matches!(system.hosted_systems[0].topology, StellarTopology::Binary(_)));
        assert!(matches!(
            system.hosted_systems[1].topology,
            StellarTopology::Single { primary_index: 0 }
        ));
        assert!(matches!(
            system.hosted_systems[2].topology,
            StellarTopology::Single { primary_index: 0 }
        ));
        assert_eq!(system.hosted_systems[1].stars[0].a, 120.0);
        assert_eq!(system.hosted_systems[2].stars[0].a, 1_200.0);
        assert_eq!(system.planetary_regions.len(), 5);
        assert_eq!(system.hosted_systems[0].planetary_regions.len(), 3);
        assert_eq!(system.hosted_systems[1].planetary_regions.len(), 1);
        assert_eq!(system.hosted_systems[2].planetary_regions.len(), 1);
        assert!(system.hosted_systems[0].planetary_regions.iter().any(|region| {
            matches!(
                region.stable_region.host,
                PlanetHost::Circumbinary {
                    hosted_system_index: 0,
                    primary_star_index: 0,
                    secondary_star_index: 1,
                }
            )
        }));
        assert!(system.hosted_systems[0].planetary_regions.iter().any(|region| {
            region.stable_region.hierarchy_host
                == Some(HierarchyHostRef::Binary {
                    node_index: 2,
                    primary_child_index: 0,
                    secondary_child_index: 1,
                })
        }));
        assert!(system.hosted_systems[1].planetary_regions.iter().all(|region| {
            matches!(
                region.stable_region.host,
                PlanetHost::Circumstellar {
                    star_ref: StarRef {
                        hosted_system_index: 1,
                        star_index: 0,
                    },
                }
            )
        }));
        assert!(system.hosted_systems[1].planetary_regions.iter().all(|region| {
            region.stable_region.hierarchy_host
                == Some(HierarchyHostRef::Star {
                    node_index: 3,
                    global_star_index: 2,
                })
        }));
        assert!(system.hosted_systems[2].planetary_regions.iter().all(|region| {
            matches!(
                region.stable_region.host,
                PlanetHost::Circumstellar {
                    star_ref: StarRef {
                        hosted_system_index: 2,
                        star_index: 0,
                    },
                }
            )
        }));
        assert!(system.hosted_systems[2].planetary_regions.iter().all(|region| {
            region.stable_region.hierarchy_host
                == Some(HierarchyHostRef::Star {
                    node_index: 5,
                    global_star_index: 3,
                })
        }));

        let hierarchy = hierarchy(&system);
        assert_eq!(hierarchy.root_index, 6);
        assert_eq!(hierarchy.nodes.len(), 7);
        assert!(matches!(
            hierarchy.nodes[0],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 0 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[1],
            StellarHierarchyNode {
                parent_index: Some(2),
                kind: StellarHierarchyNodeKind::Star { star_index: 1 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[2],
            StellarHierarchyNode {
                parent_index: Some(4),
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 0,
                    secondary_child_index: 1,
                    semi_major_axis_au: 10.0,
                    eccentricity: 0.1,
                }),
            }
        ));
        assert!(matches!(
            hierarchy.nodes[3],
            StellarHierarchyNode {
                parent_index: Some(4),
                kind: StellarHierarchyNodeKind::Star { star_index: 2 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[4],
            StellarHierarchyNode {
                parent_index: Some(6),
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 2,
                    secondary_child_index: 3,
                    semi_major_axis_au: 120.0,
                    eccentricity: 0.05,
                }),
            }
        ));
        assert!(matches!(
            hierarchy.nodes[5],
            StellarHierarchyNode {
                parent_index: Some(6),
                kind: StellarHierarchyNodeKind::Star { star_index: 3 },
            }
        ));
        assert!(matches!(
            hierarchy.nodes[6],
            StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Binary(StellarHierarchyBinary {
                    primary_child_index: 4,
                    secondary_child_index: 5,
                    semi_major_axis_au: 1_200.0,
                    eccentricity: 0.02,
                }),
            }
        ));
    }

    #[test]
    fn ambiguous_four_star_generation_returns_unresolved_stars_only_fallback() {
        let mut primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let mut secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 10.0);
        let mut tertiary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 4, 30.0);
        let mut quaternary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 6, 60.0);
        StarSystem::configure_star_orbit(&mut primary, 0.0, 0.0);
        StarSystem::configure_star_orbit(&mut secondary, 10.0, 0.6);
        StarSystem::configure_star_orbit(&mut tertiary, 30.0, 0.6);
        StarSystem::configure_star_orbit(&mut quaternary, 60.0, 0.0);

        let system = StarSystem::generate_higher_multiplicity_star_system(
            0,
            vec![primary, secondary, tertiary, quaternary],
            GenerationMode::Aggregation,
        );

        assert_eq!(system.hosted_systems.len(), 1);
        assert!(matches!(
            system.hosted_systems[0].topology,
            StellarTopology::HigherMultiplicity {
                primary_index: 0,
                star_count: 4,
            }
        ));
        assert!(system.planetary_regions.is_empty());
        assert!(system.hosted_systems[0].planetary_regions.is_empty());
        assert_eq!(system.hierarchy, Some(StellarHierarchy::unresolved(vec![0, 1, 2, 3])));
    }

    #[test]
    fn hierarchical_triple_circumbinary_region_is_limited_by_tertiary() {
        let primary = Star::new(SpectralClass::G, LuminosityClass::MainSequence, 2, 0.0);
        let secondary = Star::new(SpectralClass::K, LuminosityClass::MainSequence, 5, 8.0);
        let tertiary = Star::new(SpectralClass::M, LuminosityClass::MainSequence, 4, 140.0);
        let inner_binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: 8.0,
            eccentricity: 0.2,
            total_mass_solar: primary.mass_in_sols + secondary.mass_in_sols,
            mass_ratio: secondary.mass_in_sols / (primary.mass_in_sols + secondary.mass_in_sols),
        };
        let outer_binary = BinaryTopology {
            primary_index: 0,
            secondary_index: 1,
            semi_major_axis_au: 140.0,
            eccentricity: 0.1,
            total_mass_solar: primary.mass_in_sols + secondary.mass_in_sols + tertiary.mass_in_sols,
            mass_ratio: tertiary.mass_in_sols / (primary.mass_in_sols + secondary.mass_in_sols + tertiary.mass_in_sols),
        };

        let region = StarSystem::build_outer_limited_circumbinary_region(
            &inner_binary,
            &outer_binary,
            0,
            &primary,
            &secondary,
            GenerationMode::Aggregation,
            None,
        )
        .expect("expected hierarchical circumbinary region");

        let default_outer_au = 50.0 * (primary.mass_in_sols + secondary.mass_in_sols).cbrt();
        assert!(region.stable_region.outer_au < default_outer_au);
        assert!(region.stable_region.inner_au < region.stable_region.outer_au);
    }
}
