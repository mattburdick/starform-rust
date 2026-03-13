// src/accretion_disk.rs

use crate::random::get_random_number;
use crate::{
    accretion_parameters::get_accretion_parameters,
    body::Body,
    condensation::{CondensationThresholds, PlanetFormationInputs, RegionThermalProfile},
    consts,
    generation::GenerationMode,
    get_log_level, log,
    types::MassType,
};
use std::{collections::VecDeque, fmt};

//---------------------------  Band  ------------------------------------------
// This represents a band of dust/gas in the accretion disk
#[derive(Debug, Clone, Copy)]
pub struct Band {
    dust_present: bool,
    gas_present: bool,
    inner_edge: f64,
    outer_edge: f64,
}

impl fmt::Display for Band {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let distance = (self.outer_edge - self.inner_edge).round() as usize;
        let numchars = if distance >= 4 { distance - 4 } else { 1 }; // Subtract 4 to account for the two characters at the ends

        // Use "-" if the band has dust and gas. Otherwise, use "." (just gas)
        let output_char = if self.dust_present { "-" } else { "'" };

        // Repeat the chosen character for each AU between the previous planet and this band
        if distance > 0 {
            write!(f, "{}{:.1}", output_char.repeat(numchars), self.outer_edge)?;
        }

        Ok(())
    }
}

impl Band {
    pub fn new(inner_edge: f64, outer_edge: f64) -> Self {
        Band {
            dust_present: true,
            gas_present: true,
            inner_edge,
            outer_edge,
        }
    }

    pub fn dust_present(&self) -> bool {
        self.dust_present
    }

    pub fn gas_present(&self) -> bool {
        self.gas_present
    }

    pub fn inner_edge(&self) -> f64 {
        self.inner_edge
    }

    pub fn outer_edge(&self) -> f64 {
        self.outer_edge
    }

    fn can_merge_with(&self, other: &Band) -> bool {
        self.dust_present == other.dust_present && self.gas_present == other.gas_present
    }
    fn insert(bands: &mut Vec<Band>, new_band: Band) {
        // Find the first position where 'inner_edge' is greater
        let index = bands
            .iter()
            .position(|x| new_band.inner_edge < x.inner_edge)
            .unwrap_or(bands.len());

        bands.insert(index, new_band);
    }

    fn merge_neighbors(bands: &mut Vec<Band>) {
        let mut i = 0;
        while i < bands.len() - 1 {
            let mut merged_count = 0;
            while i + 1 < bands.len() {
                if bands[i].can_merge_with(&bands[i + 1]) {
                    bands[i].outer_edge = bands[i + 1].outer_edge; // Assume other's outer_edge is always greater
                    bands.remove(i + 1); // Remove the merged band
                    merged_count += 1;
                } else {
                    // Found the next outer neighbor that can't be merged with, so stop trying to merge
                    break;
                }
            }
            i += merged_count + 1;
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct EmbryoCell {
    orbital_radius_au: f64,
    orbital_eccentricity: f64,
    inner_edge_au: f64,
    outer_edge_au: f64,
}

//---------------------------  AccretionStepResult  ---------------------------

/// Record of a single collision that occurred during coalescence.
#[derive(Debug, Clone)]
pub struct CollisionRecord {
    /// Semi-major axis of the accreting body before the collision.
    pub body_pre_a_au: f64,
    /// Whether the accreting body was a gas giant before the collision.
    pub body_was_gas_giant: bool,
    /// Radius of the accreting body in km before the collision.
    pub body_radius_in_km: f64,
    /// Semi-major axis of the absorbed neighbor before the collision.
    pub neighbor_pre_a_au: f64,
    /// Whether the absorbed neighbor was a gas giant.
    pub neighbor_was_gas_giant: bool,
    /// Radius of the absorbed neighbor in km (for rendering).
    pub neighbor_radius_in_km: f64,
    /// Semi-major axis of the merged body after the collision.
    pub result_a_au: f64,
}

/// The result of a single aggregation step (one protoplanet injection attempt).
#[derive(Debug, Clone)]
pub struct AccretionStepResult {
    pub injection_radius_au: f64,
    pub injection_eccentricity: f64,
    pub outcome: StepOutcome,
    /// Collisions that occurred during coalescence (empty if none).
    pub collisions: Vec<CollisionRecord>,
}

/// Whether a protoplanet injection formed a body or failed.
#[derive(Debug, Clone, PartialEq)]
pub enum StepOutcome {
    BodyFormed,
    FailedInjection,
}

//---------------------------  AccretionDisk  ---------------------------------

#[derive(Debug, Clone)]
pub struct AccretionDisk {
    pub central_mass_in_sols: f64,
    pub luminosity_in_sols: f64, // The luminosity of the nearby star in solar luminosities
    pub thermal_profile: RegionThermalProfile,
    pub planet_inner_bound: f64, // Inner limit at which a body can exist in orbit about the central mass
    pub planet_outer_bound: f64, // Outer limit at which a body can exist in orbit about the central mass
    pub disk_inner_bound: f64,   // Inner limit of the accretion disk
    pub disk_outer_bound: f64,   // Outer limit of the accretion disk
    pub bands: Vec<Band>,        // The bands of dust and gas comprising the accretion disk
    pub bodies: Vec<Body>,       // The bodies in the accretion disk
    pub dust_left: bool,         // If true, dust remains to be accreted
    /// Cached gas-to-dust ratio, read once from ACCRETION_PARAMETERS at construction.
    /// Avoids repeated mutex locks in the hot collect_dust inner loop.
    cached_gas_to_dust_ratio: f64,
    /// Cached dust density coefficient from ACCRETION_PARAMETERS.
    cached_dust_density_coeff: f64,
}

impl fmt::Display for AccretionDisk {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut previous_a = 0.0;
        for body in &self.bodies {
            let distance = (body.a - previous_a).round() as usize;
            previous_a = body.a;
            if distance > 0 {
                write!(
                    f,
                    "{}{:.1}",
                    " ".repeat(distance),
                    body.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES
                )?;
            }
        }
        writeln!(f)?;

        for body in &self.bodies {
            let distance = body.a.round() as usize;
            let body_char = match body.mass_type {
                MassType::Moon => "o",
                MassType::Planet => "o",
                MassType::GasGiant => "0",
                MassType::Star => "*",
            };
            if distance > 0 {
                writeln!(f, "{}{body_char}", " ".repeat(distance))?;
            }
        }
        writeln!(f)?;

        for band in &self.bands {
            write!(f, "{}", band)?;
        }
        writeln!(f)?;

        Ok(())
    }
}

impl Default for AccretionDisk {
    fn default() -> Self {
        Self::new(&Body::default(), 0.0)
    }
}

impl AccretionDisk {
    const SEMI_ANALYTIC_BASE_SPACING_RATIO: f64 = 1.65;
    const SEMI_ANALYTIC_SPACING_JITTER: f64 = 0.35;
    const SEMI_ANALYTIC_MIN_SPACING_RATIO: f64 = 1.18;
    const SEMI_ANALYTIC_MAX_SPACING_RATIO: f64 = 2.60;
    const SEMI_ANALYTIC_ECCENTRICITY_SCALE: f64 = 0.35;
    const SEMI_ANALYTIC_ISOLATION_SCALE: f64 = 6.0e3;
    const SEMI_ANALYTIC_MIN_SOLID_FRACTION: f64 = 0.02;

    pub fn new_with_planet_bounds(
        central_mass_in_sols: f64,
        luminosity_in_sols: f64,
        central_mass_type: MassType,
        planet_inner_bound: f64,
        planet_outer_bound: f64,
    ) -> Self {
        let clamped_inner_bound = planet_inner_bound.max(0.0);
        let clamped_outer_bound = planet_outer_bound.max(0.0);
        let (disk_inner_bound, _, _, _) =
            Body::gravitational_effect_limits(clamped_inner_bound, 0.0, consts::PROTOPLANET_MASS);
        let (_, raw_disk_outer_bound, _, _) =
            Body::gravitational_effect_limits(clamped_outer_bound, 0.0, consts::PROTOPLANET_MASS);

        let disk_outer_bound =
            raw_disk_outer_bound.min(Self::stell_dust_limit(central_mass_in_sols, 0.0, central_mass_type));

        let has_viable_band = clamped_outer_bound > clamped_inner_bound && disk_outer_bound > disk_inner_bound;
        let mut bands: Vec<Band> = Vec::new();
        if has_viable_band {
            Band::insert(&mut bands, Band::new(disk_inner_bound, disk_outer_bound));
        }

        let (cached_gas_to_dust_ratio, cached_dust_density_coeff) = {
            let params = crate::accretion_parameters::ACCRETION_PARAMETERS.lock().unwrap();
            (
                ((1.0 / params.percent_dust_in_cloud) * 100.0).min(1000.0),
                params.dust_density_coefficient,
            )
        };

        AccretionDisk {
            luminosity_in_sols,
            central_mass_in_sols,
            thermal_profile: RegionThermalProfile::from_host_luminosity(luminosity_in_sols),
            planet_inner_bound: clamped_inner_bound,
            planet_outer_bound: clamped_outer_bound,
            disk_inner_bound,
            disk_outer_bound,
            bands,
            bodies: Vec::new(),
            dust_left: has_viable_band,
            cached_gas_to_dust_ratio,
            cached_dust_density_coeff,
        }
    }

    /// Generates a random eccentricity value for an orbital body.
    ///
    /// This function simulates a distribution of eccentricities where values are more likely to be closer to 1.
    /// It generates a random floating-point number between 0.0001 and 1.0, inclusive, to ensure non-zero eccentricity.
    /// The resulting eccentricity is adjusted using a power function to skew towards higher values, simulating a more
    /// realistic distribution of orbital eccentricities found in nature.
    ///
    /// The constant `ECCENTRICITY_COEFF` is used to control the skewness of the distribution towards higher values.
    /// This constant should ideally reflect the typical distribution seen in celestial bodies.
    ///
    /// # Returns
    /// A `f64` representing the eccentricity, adjusted to be typically closer to 1.
    ///
    /// # Examples
    /// ```rust
    /// use starform_rust::accretion_disk::AccretionDisk;
    ///
    /// let eccentricity = AccretionDisk::random_eccentricity();
    /// println!("Generated eccentricity: {}", eccentricity);
    /// ```
    pub fn random_eccentricity() -> f64 {
        let random_number_between_0_and_1: f64 = get_random_number(0.0001..=1.0);
        1.0 - (1.0 - random_number_between_0_and_1).powf(consts::ECCENTRICITY_COEFF)
    }

    // Calculates the outer limit of gas & dust accretion disks from a planet or star
    fn stell_dust_limit(central_mass_in_sols: f64, dist_from_primary: f64, central_mass_type: MassType) -> f64 {
        let mut outer_limit = 200.0 * central_mass_in_sols.cbrt();
        if central_mass_type == MassType::Star {
            return outer_limit;
        }

        // Otherwise the central mass is a planet
        outer_limit /= 125.0;
        let primary_effect = dist_from_primary.powi(2);
        if primary_effect <= 1.0 {
            outer_limit *= primary_effect;
        }
        outer_limit
    }

    #[allow(dead_code)]
    fn dust_available(&self, inside_range: f64, outside_range: f64) -> bool {
        for band in &self.bands {
            if !band.dust_present {
                continue;
            }
            if (outside_range.min(band.outer_edge) - inside_range.max(band.inner_edge)) > 0.0 {
                return true;
            }
        }

        false
    }

    /// Calculates the volume of a cylindrical shell (annular cylinder) given its dimensions.
    ///
    /// The volume is determined based on the inner and outer radii of the cylinder and the sum of xa and xp,
    /// the gravitational influence of a mass at aphelion and perihelion respectively. All units are in AU.
    ///
    /// # Parameters
    /// - `r_inner`: The inner radius of the cylindrical shell in units (e.g., meters).
    /// - `r_outer`: The outer radius of the cylindrical shell in the same units as `r_inner`.
    /// - `xp`: Radius of a mass's gravitational influence at perihelion.
    /// - `xa`: Radius of a mass's gravitational influence at aphelion.
    ///
    /// # Returns
    /// Returns the volume of the cylindrical shell in cubic units, calculated using the formula:
    /// \[
    /// V=\pi \cdot (x_{p} + x_{a}) \cdot (r_{a}^2 - r_{p}^2)
    /// \\
    /// \text{Where a=aphelion and p=perihelion}
    /// \]
    ///
    /// # Examples
    /// ```rust
    /// use starform_rust::accretion_disk::AccretionDisk;
    ///
    /// let volume = AccretionDisk::volume(2.0, 3.0, 1.0, 1.0);
    /// println!("Volume of the cylindrical shell: {}", volume);
    /// ```
    ///
    /// This calculates the volume of a cylindrical shell with an inner radius of 2 units, an outer radius of 3 units,
    /// and a total height of 2 units.
    pub fn volume(r_inner: f64, r_outer: f64, xp: f64, xa: f64) -> f64 {
        // Total height of the cylindrical shell
        let height = xa + xp;

        // Calculating the volume
        std::f32::consts::PI as f64 * height * (r_outer.powi(2) - r_inner.powi(2))
    }

    /// Determines if a new body is dynamically overpacked with any existing
    /// body using the same orbital-ratio and mutual-Hill-sphere criteria
    /// applied by the post-accretion spacing cleanup.  When two bodies
    /// violate either criterion they are merged during accretion, preserving
    /// mass in the system instead of discarding the smaller body later.
    ///
    /// # Returns
    /// `(true, index)` when the closest overpacked neighbor is at `index`,
    /// or `(false, 0)` when no collision is detected.
    fn find_collision(&self, new_body: &Body) -> (bool, usize) {
        let mut closest_neighbor = 0;
        let mut found_collision = false;
        let mut closest_approach = f64::INFINITY;

        for index in 0..self.bodies.len() {
            let body = &self.bodies[index];

            if !Self::bodies_are_overpacked(self.central_mass_in_sols, body, new_body) {
                continue;
            }

            let separation = (body.a - new_body.a).abs();
            if !found_collision || separation < closest_approach {
                closest_neighbor = index;
                closest_approach = separation;
                found_collision = true;
            }
        }

        (found_collision, closest_neighbor)
    }

    /// Two bodies are "overpacked" if their orbital semi-major axis ratio is
    /// below `MIN_PLANETARY_LOG_SPACING_RATIO` *or* they are separated by
    /// fewer than `MIN_PLANETARY_MUTUAL_HILL_SPACING` mutual Hill radii.
    fn bodies_are_overpacked(host_mass_solar: f64, a: &Body, b: &Body) -> bool {
        let (inner, outer) = if a.a <= b.a { (a, b) } else { (b, a) };

        let ratio = if inner.a > 0.0 {
            outer.a / inner.a
        } else {
            f64::INFINITY
        };
        if ratio < consts::MIN_PLANETARY_LOG_SPACING_RATIO {
            return true;
        }

        let combined_mass_solar = (inner.mass_in_sols + outer.mass_in_sols).max(0.0);
        if host_mass_solar <= 0.0 || combined_mass_solar <= 0.0 {
            return false;
        }

        let average_a = (inner.a + outer.a) * 0.5;
        let mutual_hill_radius = average_a * (combined_mass_solar / (3.0 * host_mass_solar)).cbrt();
        let required_spacing = consts::MIN_PLANETARY_MUTUAL_HILL_SPACING * mutual_hill_radius;
        let actual_spacing = (outer.a - inner.a).max(0.0);

        actual_spacing < required_spacing
    }

    fn refresh_body_accretion_state(&self, body: &mut Body) {
        body.local_dust_density = self.cached_dust_density_coeff
            * self.central_mass_in_sols.sqrt()
            * f64::exp(-consts::ALPHA * body.a.cbrt());
        body.critical_mass_limit = body.critical_limit(self.luminosity_in_sols);
    }

    fn sample_planet_formation_inputs(&self, orbital_radius_au: f64) -> PlanetFormationInputs {
        self.thermal_profile
            .sample_planet_formation_inputs(orbital_radius_au, CondensationThresholds::default())
    }

    fn clear_compatibility_bands(&mut self) {
        self.bands.clear();
    }

    fn reset_generation_state(&mut self) {
        self.bodies.clear();
        self.bands.clear();

        if self.disk_outer_bound > self.disk_inner_bound && self.planet_outer_bound > self.planet_inner_bound {
            Band::insert(&mut self.bands, Band::new(self.disk_inner_bound, self.disk_outer_bound));
            self.dust_left = true;
        } else {
            self.dust_left = false;
        }
    }

    fn sample_log_uniform(min_value: f64, max_value: f64) -> f64 {
        let min_value = min_value.max(1.0e-4);
        let max_value = max_value.max(min_value);
        if (max_value - min_value).abs() <= f64::EPSILON {
            return min_value;
        }

        let min_log = min_value.log10();
        let max_log = max_value.log10();
        10_f64.powf(get_random_number(min_log..=max_log))
    }

    fn sample_embryo_step_ratio() -> f64 {
        let jitter =
            get_random_number((1.0 - Self::SEMI_ANALYTIC_SPACING_JITTER)..=(1.0 + Self::SEMI_ANALYTIC_SPACING_JITTER));
        (Self::SEMI_ANALYTIC_BASE_SPACING_RATIO * jitter).clamp(
            Self::SEMI_ANALYTIC_MIN_SPACING_RATIO,
            Self::SEMI_ANALYTIC_MAX_SPACING_RATIO,
        )
    }

    fn generate_log_spaced_embryo_cells(&self) -> Vec<EmbryoCell> {
        if self.planet_outer_bound <= self.planet_inner_bound {
            return Vec::new();
        }

        let mut orbital_radii_au = Vec::new();
        let initial_max = (self.planet_inner_bound.max(0.01) * Self::SEMI_ANALYTIC_BASE_SPACING_RATIO.sqrt())
            .min(self.planet_outer_bound)
            .max(self.planet_inner_bound.max(0.01));
        let mut orbital_radius_au = Self::sample_log_uniform(self.planet_inner_bound.max(0.01), initial_max);

        while orbital_radius_au <= self.planet_outer_bound {
            orbital_radii_au.push(orbital_radius_au);
            let step_ratio = Self::sample_embryo_step_ratio();
            let next_orbital_radius_au = orbital_radius_au * step_ratio;
            if next_orbital_radius_au <= orbital_radius_au {
                break;
            }
            orbital_radius_au = next_orbital_radius_au;
        }

        if orbital_radii_au.is_empty() {
            orbital_radii_au.push((self.planet_inner_bound + self.planet_outer_bound) / 2.0);
        }

        orbital_radii_au
            .iter()
            .enumerate()
            .filter_map(|(index, orbital_radius_au)| {
                let inner_edge_au = if index == 0 {
                    self.planet_inner_bound
                } else {
                    (orbital_radii_au[index - 1] * orbital_radius_au)
                        .sqrt()
                        .max(self.planet_inner_bound)
                };
                let outer_edge_au = if index + 1 == orbital_radii_au.len() {
                    self.planet_outer_bound
                } else {
                    (orbital_radius_au * orbital_radii_au[index + 1])
                        .sqrt()
                        .min(self.planet_outer_bound)
                };

                if outer_edge_au <= inner_edge_au {
                    return None;
                }

                Some(EmbryoCell {
                    orbital_radius_au: Self::sample_log_uniform(inner_edge_au, outer_edge_au),
                    orbital_eccentricity: (Self::random_eccentricity() * Self::SEMI_ANALYTIC_ECCENTRICITY_SCALE)
                        .min(0.35),
                    inner_edge_au,
                    outer_edge_au,
                })
            })
            .collect()
    }

    fn estimate_local_solid_surface_density(&self, body: &Body, formation_inputs: PlanetFormationInputs) -> f64 {
        let solids_fraction = (1.0 - formation_inputs.condensation_fractions.gas).clamp(0.0, 1.0);
        if solids_fraction < Self::SEMI_ANALYTIC_MIN_SOLID_FRACTION {
            return 0.0;
        }

        body.local_dust_density * solids_fraction
    }

    fn estimate_isolation_mass(&self, cell: EmbryoCell, formation_inputs: PlanetFormationInputs, body: &Body) -> f64 {
        let solid_surface_density = self.estimate_local_solid_surface_density(body, formation_inputs);
        if solid_surface_density <= 0.0 {
            return 0.0;
        }

        let annulus_area = (cell.outer_edge_au.powi(2) - cell.inner_edge_au.powi(2)).max(0.0);
        if annulus_area <= 0.0 {
            return 0.0;
        }

        let icy_enhancement = 1.0
            + formation_inputs.condensation_fractions.water_ice
            + formation_inputs.condensation_fractions.volatile_ices;
        let orbital_scale = cell.orbital_radius_au.max(0.1).sqrt();
        let host_mass_scale = self.central_mass_in_sols.max(0.1).powf(0.15);

        Self::SEMI_ANALYTIC_ISOLATION_SCALE * solid_surface_density * annulus_area * orbital_scale * icy_enhancement
            / host_mass_scale
    }

    fn estimate_envelope_mass(
        &self,
        core_mass_in_sols: f64,
        critical_mass_limit: f64,
        formation_inputs: PlanetFormationInputs,
    ) -> f64 {
        if critical_mass_limit <= 0.0 || core_mass_in_sols <= 0.0 {
            return 0.0;
        }

        let gas_fraction = formation_inputs.condensation_fractions.gas;
        if gas_fraction <= 0.0 {
            return 0.0;
        }

        let accretion_parameters = get_accretion_parameters();
        let gas_to_dust_ratio =
            ((100.0 / accretion_parameters.percent_dust_in_cloud.max(0.1)) - 1.0).clamp(0.0, 1000.0);
        let giant_core_eligibility = ((core_mass_in_sols / critical_mass_limit) - 0.75).clamp(0.0, 1.5);

        core_mass_in_sols * giant_core_eligibility * gas_fraction * ((gas_to_dust_ratio + 1.0) / 50.0).sqrt() * 0.65
    }

    fn finalize_generated_body(&self, body: &mut Body, force_gas_giant: bool) {
        let formation_inputs = self.sample_planet_formation_inputs(body.a);
        let gas_available = formation_inputs.condensation_fractions.gas > 0.05;
        body.mass_type = if force_gas_giant
            || (gas_available && body.critical_mass_limit > 0.0 && body.mass_in_sols >= body.critical_mass_limit)
        {
            MassType::GasGiant
        } else {
            MassType::Planet
        };
        body.initialize_from_formation_inputs(formation_inputs);
    }

    fn estimate_candidate_body(&self, cell: EmbryoCell) -> Option<Body> {
        let mut body = Body::new(
            cell.orbital_radius_au,
            cell.orbital_eccentricity,
            consts::PROTOPLANET_MASS,
            MassType::Planet,
            0.0,
            self.central_mass_in_sols,
            self.luminosity_in_sols,
            None,
        );
        let formation_inputs = self.sample_planet_formation_inputs(cell.orbital_radius_au);
        let core_mass_in_sols = self.estimate_isolation_mass(cell, formation_inputs, &body);
        if core_mass_in_sols <= consts::TRIVIAL_MASS {
            return None;
        }

        body.mass_in_sols = core_mass_in_sols;
        self.refresh_body_accretion_state(&mut body);
        let envelope_mass_in_sols =
            self.estimate_envelope_mass(core_mass_in_sols, body.critical_mass_limit, formation_inputs);
        body.mass_in_sols += envelope_mass_in_sols;
        self.refresh_body_accretion_state(&mut body);
        self.finalize_generated_body(&mut body, false);
        Some(body)
    }

    fn coalesce_generated_body(&mut self, mut body: Body) -> Body {
        let mut force_gas_giant = body.mass_type == MassType::GasGiant;
        loop {
            let (found_collision, closest_neighbor) = self.find_collision(&body);
            if !found_collision {
                break;
            }

            let neighbor = self.bodies.remove(closest_neighbor);
            force_gas_giant = force_gas_giant || neighbor.mass_type == MassType::GasGiant;

            body.collide(&neighbor);
            self.refresh_body_accretion_state(&mut body);
        }

        self.finalize_generated_body(&mut body, force_gas_giant);
        body
    }

    #[allow(dead_code)]
    fn coalesce_body(&mut self, mut body: Body) -> (Body, Vec<CollisionRecord>) {
        let mut collisions = Vec::new();
        loop {
            let (found_collision, closest_neighbor) = self.find_collision(&body);
            if !found_collision {
                break;
            }

            let body_pre_a = body.a;
            let body_was_gas_giant = body.mass_type == MassType::GasGiant;
            let body_radius_in_km = body.radius_in_km;
            let neighbor = self.bodies.remove(closest_neighbor);
            let neighbor_pre_a = neighbor.a;
            let neighbor_was_gas_giant = neighbor.mass_type == MassType::GasGiant;
            let neighbor_radius_in_km = neighbor.radius_in_km;

            let should_remain_gas_giant = body.mass_type == MassType::GasGiant || neighbor_was_gas_giant;

            body.collide(&neighbor);
            if should_remain_gas_giant {
                body.mass_type = MassType::GasGiant;
            }

            collisions.push(CollisionRecord {
                body_pre_a_au: body_pre_a,
                body_was_gas_giant,
                body_radius_in_km,
                neighbor_pre_a_au: neighbor_pre_a,
                neighbor_was_gas_giant,
                neighbor_radius_in_km,
                result_a_au: body.a,
            });

            self.refresh_body_accretion_state(&mut body);
            body = self.accrete_dust(body);
        }

        if body.critical_mass_limit > 0.0 && body.mass_in_sols >= body.critical_mass_limit {
            body.mass_type = MassType::GasGiant;
        }

        let formation_inputs = self.sample_planet_formation_inputs(body.a);
        body.initialize_from_formation_inputs(formation_inputs);
        (body, collisions)
    }

    /// Simulates the process of dust and gas accretion for a celestial body in an accretion disk.
    ///
    /// This function iteratively simulates the process of a celestial body, represented by the `Body` struct,
    /// sweeping up dust and gas from its orbital path. It repeatedly calls the `collect_dust` method, which
    /// increases the body's mass based on the amount of dust and gas it can accrete during each iteration.
    /// The process continues until the additional mass gained in an iteration is less than 0.01% of the total mass,
    /// indicating that accretion has slowed down significantly.
    ///
    /// During the accretion process, if the body's mass exceeds the `critical_mass_limit`, it is classified as a
    /// `GasGiant`. After the loop concludes, the function checks the entire accretion disk to determine if any
    /// dust is left, updating the `dust_left` flag accordingly.
    ///
    /// Detailed logging is provided throughout the accretion process, including the final mass, orbital distance,
    /// number of iterations taken, and gravitational effect limits if the logging level is sufficiently high.
    ///
    /// # Arguments
    /// * `body` - A mutable reference to the `Body` being accreted upon.
    ///
    /// # Returns
    /// Returns the updated `Body` struct after accretion.
    ///
    /// # Examples
    /// ```rust,ignore
    /// let mut body = Body {
    ///     mass_in_sols: 0.05,
    ///     critical_mass_limit: 0.08,
    ///     ... // other fields
    /// };
    /// let accreted_body = accretion_disk.accrete_dust(body);
    /// println!("New mass in solar masses: {}", accreted_body.mass_in_sols);
    /// ```
    #[allow(dead_code)]
    fn accrete_dust(&mut self, mut body: Body) -> Body {
        let mut loop_count = 0;
        loop {
            loop_count += 1;
            let old_mass = body.mass_in_sols;
            body.mass_in_sols = self.collect_dust(&body);
            if body.mass_in_sols >= body.critical_mass_limit {
                body.mass_type = MassType::GasGiant;
            }

            if (body.mass_in_sols - old_mass) <= (0.0001 * old_mass) {
                break;
            }
        }

        // Check if there is any dust remaining
        self.dust_left = false;
        for band in &self.bands {
            if band.dust_present {
                self.dust_left = true;
                break;
            }
        }

        let log_level = *get_log_level!();
        log!(
            log_level,
            1,
            "Accreted a {} of {:.3} Earth masses at {:.3} AU over {loop_count} iterations (critical mass limit: {:.3} Earth masses)",
            body.mass_type,
            body.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES,
            body.a,
            body.critical_mass_limit * consts::SUN_MASS_IN_EARTH_MASSES
        );

        if log_level >= 2 {
            let (r_inner, r_outer, xp, xa) = Body::gravitational_effect_limits(body.a, body.e, body.mass_in_sols);
            log!(
                log_level,
                2,
                "r_inner={:.2}, r_outer={:.2}, xp={:.2}, xa={:.2}",
                r_inner,
                r_outer,
                xp,
                xa
            );
        }

        body
    }

    /// Simulates the accretion process by a celestial body collecting dust and gas within its gravitational influence.
    ///
    /// This function iterates over each dust and gas band in the accretion disk, determining if the band falls within
    /// the gravitational effect limits (from 'r_inner' to 'r_outer') of the provided celestial body (`protoplanet`). If a band
    /// is within the influence range, it checks if the protoplanet's mass is sufficient to collect dust and, if the mass
    /// exceeds the critical mass, gas as well. Bands within the range but without sufficient mass or lacking dust/gas
    /// are skipped.
    ///
    /// The function updates the mass of the protoplanet by adding the mass of the dust and gas collected. It also
    /// manages the creation of new bands for gaps created by partial band accretion. These new bands are stored
    /// temporarily and merged into the main list of bands at the end of the operation.
    ///
    /// # Parameters
    /// * `protoplanet`: A reference to the `Body` struct representing the accreting celestial body.
    ///
    /// # Returns
    /// Returns the new total mass of the protoplanet after the accretion process.
    ///
    /// # Side Effects
    /// * Modifies the list of bands within the `AccretionDisk` to reflect dust and gas depletion and new band formation.
    /// * Updates the `dust_left` flag of the `AccretionDisk` based on remaining dust in any band.
    ///
    /// # Example
    /// ```rust,ignore
    /// let mut accretion_disk = AccretionDisk::new(...);
    /// let mut protoplanet = Body {
    ///     mass_in_sols: 0.05,
    ///     a: 5.0,
    ///     e: 0.1,
    ///     ... // other fields
    /// };
    /// let new_mass = accretion_disk.collect_dust(&protoplanet);
    /// println!("New mass after accretion: {}", new_mass);
    /// ```
    pub fn collect_dust(&mut self, protoplanet: &Body) -> f64 {
        // The injected mass has inner and outer effect limits based on size of the mass, its orbit, and eccentricity.
        // We will accumulate mass from all the dust and gas bands into "mass".

        // Calculate the inner and outer effect limits of the injected mass
        // The inner and outer effect limits are the distances from the primary star at which the mass can affect the dust bands.
        // xa and xp represent the aphelion and parahelion distances from the orbital plane at which the mass can affect the dust bands.
        let (r_inner, r_outer, xp, xa) =
            Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass_in_sols);

        let mut new_mass = protoplanet.mass_in_sols;

        // Sweeping out dust from part of a band may create gaps represented as new bands that need to be added to the list.
        // We'll store them in a queue and add them after we've iterated over all the bands.
        let mut new_bands: VecDeque<Band> = VecDeque::new();
        for band in &mut self.bands {
            let outer_gap = band.outer_edge - r_outer;
            let inner_gap = band.inner_edge - r_inner;

            if r_outer <= band.inner_edge || r_inner >= band.outer_edge {
                // The bandwidth swept out by the mass is either inside the inner edge of this band or outside the outer edge.
                // Skip this band.
                continue;
            }

            if (!protoplanet.collects_gas(new_mass) && !band.dust_present)
                || (protoplanet.collects_gas(new_mass) && !band.gas_present && !band.dust_present)
            {
                // The mass is too small to pick up gas and the band has no dust.
                continue;
            }

            // Start by assuming the entire bandwidth is swept up. We'll subtract the gaps later, recalling that they can be negative.
            if inner_gap < 0.0 {
                // The effect limit is outside the inner edge of the band. Reduce the swept bandwidth by the inner gap and create
                // a new band representing the gap.

                // Create a new band for the gap.
                let mut gap_band = Band::new(band.inner_edge, r_inner);
                gap_band.dust_present = band.dust_present;
                gap_band.gas_present = band.gas_present;
                new_bands.push_back(gap_band);

                // Move the inner edge out
                band.inner_edge = r_inner;
            }
            if outer_gap > 0.0 {
                // Subtract the outer gap from the swept bandwidth and create a new band for the gap.
                let mut gap_band = Band::new(r_outer, band.outer_edge);
                gap_band.dust_present = band.dust_present;
                gap_band.gas_present = band.gas_present;
                new_bands.push_back(gap_band);

                // Move the outer edge in
                band.outer_edge = r_outer;
            }

            band.dust_present = false;
            band.gas_present = if protoplanet.collects_gas(new_mass) {
                false
            } else {
                band.gas_present
            };

            let volume = Self::volume(r_inner, r_outer, xp, xa);
            let ratio = self.cached_gas_to_dust_ratio;
            let density = ratio * protoplanet.local_dust_density
                / (1.0 + (protoplanet.critical_mass_limit / new_mass).sqrt() * (ratio - 1.0));
            new_mass += volume * density;
        }

        // Apply all changes here
        for band in new_bands {
            Band::insert(&mut self.bands, band);
        }

        Band::merge_neighbors(&mut self.bands);

        // TODO: do we neeed to update the bands in the body?
        // // Now update `dust_left` based on whether any band still has dust.
        // self.dust_left = self.bands.iter().any(|band| band.dust_present);

        new_mass
    }

    /// Creates a new `AccretionDisk` around a primary body with specified luminosity.
    ///
    /// This method initializes an `AccretionDisk` with bounds calculated based on the gravitational
    /// and luminous characteristics of the primary body. It considers both the physical limits imposed
    /// by the body's Roche limit and the effects of its luminosity on the surrounding dust and gas.
    ///
    /// # Parameters:
    /// - `primary`: Reference to the central `Body` object around which the accretion disk forms.
    /// - `luminosity_in_sols`: The luminosity of the primary body in solar luminosities.
    ///
    /// # Returns:
    /// - `Self`: An instance of `AccretionDisk` with all fields initialized, including calculated
    ///   inner and outer bounds for planets and the dust/gas disk.
    ///
    /// # Calculations:
    /// - `planet_inner_bound`: Calculated as the minimum of the Roche limit and a factor of the primary's
    ///   mass if the body has significant luminosity.
    /// - `planet_outer_bound`: Set to a function of the primary's mass, typically 50 times the cube root
    ///   of the mass.
    /// - `disk_inner_bound` and `disk_outer_bound`: Define the range of the dust and gas disk. The inner
    ///   bound is set by the gravitational effect limits at the inner planet bound with zero eccentricity.
    ///   The outer bound is similarly set but adjusted based on stellar properties and might be truncated
    ///   based on the type of star.
    ///
    /// # Example Usage:
    /// ```rust,ignore
    /// let primary_star = Body {
    ///     mass_in_sols: 1.0,
    ///     a: 0.0,
    ///     e: 0.0,
    ///     ...
    /// };
    /// let accretion_disk = AccretionDisk::new(&primary_star, 1.0);
    /// ```
    ///
    /// The created accretion disk includes a vector of bands representing the initial distribution of dust
    /// and gas based on the calculated boundaries.
    pub fn new(primary: &Body, luminosity_in_sols: f64) -> Self {
        // Choose the lower of the Roche limit for a fluid satellite and a limit due to the luminosity of the central star (if any)
        let mut planet_inner_bound = primary.roche_limit_in_au();
        if luminosity_in_sols > 0.0 {
            planet_inner_bound = planet_inner_bound.min(0.3 * primary.mass_in_sols.cbrt());
        }

        let planet_outer_bound = 50.0 * primary.mass_in_sols.cbrt();

        Self::new_with_planet_bounds(
            primary.mass_in_sols,
            luminosity_in_sols,
            primary.mass_type,
            planet_inner_bound,
            planet_outer_bound,
        )
    }

    /// Simulates the accretion process within an accretion disk by creating and evolving protoplanets.
    ///
    /// This function repeatedly checks for dust presence in the accretion disk and initiates the
    /// accretion process by creating protoplanets in bands where dust is available. Each protoplanet
    /// undergoes dust and gas accumulation based on its gravitational influence range. The function
    /// checks for possible collisions with existing bodies in the system and handles them accordingly.
    ///
    /// # Workflow:
    /// 1. Determine if there's any dust left to accrete.
    /// 2. Randomly determine the eccentricity for the protoplanet.
    /// 3. Locate a dust band that still contains dust.
    /// 4. Randomly choose an orbital distance within the dust band's effective gravitational range.
    /// 5. Instantiate a protoplanet at the chosen location with initial properties.
    /// 6. Calculate the protoplanet's gravitational effect limits and ensure there's dust available in that range.
    /// 7. Accrete dust and gas onto the protoplanet until the mass increment becomes trivial.
    /// 8. If the protoplanet collides with another body, merge them and reaccrete dust as needed.
    /// 9. If no collision occurs and the protoplanet achieves or exceeds the critical mass limit, classify it as a Gas Giant.
    /// 10. Add the protoplanet to the system if it remains after potential collisions.
    ///
    /// # Returns
    /// Returns a mutable reference to the `AccretionDisk` to allow for chaining and further modifications.
    ///
    /// # Example
    /// ```rust,ignore
    /// let mut accretion_disk = AccretionDisk::new(...);
    /// accretion_disk.accrete();
    /// println!("Accretion process completed: {}", accretion_disk);
    /// ```
    ///
    /// # Panics
    /// Panics if no dust band with dust is found, indicating an error in managing the state of the disk.
    /// Also, it will panic if the calculated bounds for the planet's orbit are out of expected range,
    /// which could indicate incorrect initialization or corruption of disk state.
    /// Performs one iteration of the aggregation accretion loop: picks a random
    /// dust band, injects a protoplanet, accretes dust/gas, coalesces collisions,
    /// and returns a step result.
    ///
    /// After this call, `self.bands`, `self.bodies`, and `self.dust_left` are
    /// updated to reflect the completed step.
    fn aggregation_step(&mut self) -> AccretionStepResult {
        let orbital_eccentricity = Self::random_eccentricity();
        let band = self
            .bands
            .iter()
            .find(|band| band.dust_present)
            .expect("AccretionDisk.aggregation_step: unable to find a band with dust");

        let bound1 = band
            .inner_edge
            .max(self.planet_inner_bound)
            .min(self.planet_outer_bound);
        let bound2 = band
            .outer_edge
            .min(self.planet_outer_bound)
            .max(self.planet_inner_bound);
        let orbital_radius_au = get_random_number(bound1..=bound2);

        let mut protoplanet = Body::new(
            orbital_radius_au,
            orbital_eccentricity,
            consts::PROTOPLANET_MASS,
            MassType::Planet,
            0.0,
            self.central_mass_in_sols,
            self.luminosity_in_sols,
            None,
        );
        let (eff_inner_bound, eff_outer_bound, _, _) =
            Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass_in_sols);

        if !self.dust_available(eff_inner_bound, eff_outer_bound) {
            panic!("ERROR: found no dust in range of protoplanet");
        }

        protoplanet = self.accrete_dust(protoplanet);
        if protoplanet.is_trivial_mass() {
            let log_level = *get_log_level!();
            log!(
                log_level,
                1,
                "Failed injection at {:.3} AU (trivial mass after accretion: {:.6} Earth masses)",
                orbital_radius_au,
                protoplanet.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES
            );
            return AccretionStepResult {
                injection_radius_au: orbital_radius_au,
                injection_eccentricity: orbital_eccentricity,
                outcome: StepOutcome::FailedInjection,
                collisions: Vec::new(),
            };
        }

        let (protoplanet, collisions) = self.coalesce_body(protoplanet);
        if protoplanet.is_trivial_mass() {
            let log_level = *get_log_level!();
            log!(
                log_level,
                1,
                "Failed injection at {:.3} AU (trivial mass after coalescence: {:.6} Earth masses)",
                orbital_radius_au,
                protoplanet.mass_in_sols * consts::SUN_MASS_IN_EARTH_MASSES
            );
            return AccretionStepResult {
                injection_radius_au: orbital_radius_au,
                injection_eccentricity: orbital_eccentricity,
                outcome: StepOutcome::FailedInjection,
                collisions,
            };
        }

        Body::insert(&mut self.bodies, protoplanet);

        AccretionStepResult {
            injection_radius_au: orbital_radius_au,
            injection_eccentricity: orbital_eccentricity,
            outcome: StepOutcome::BodyFormed,
            collisions,
        }
    }

    fn accrete_aggregation(&mut self) -> &mut Self {
        self.reset_generation_state();
        if self.planet_inner_bound > self.planet_outer_bound {
            panic!("ERROR: planet inner bound is greater than outer bound\n");
        }

        while self.dust_left {
            self.aggregation_step();
        }

        let log_level = *get_log_level!();
        log!(
            log_level,
            1,
            "Generated {} aggregated bodies between {:.2} AU and {:.2} AU",
            self.bodies.len(),
            self.planet_inner_bound,
            self.planet_outer_bound
        );
        if log_level >= 2 && !self.bodies.is_empty() {
            log!(log_level, 2, "{}", self);
        }

        self
    }

    fn accrete_semi_analytic(&mut self) -> &mut Self {
        self.bodies.clear();
        self.clear_compatibility_bands();

        if self.planet_inner_bound > self.planet_outer_bound {
            panic!("ERROR: planet inner bound is greater than outer bound\n");
        }

        for embryo_cell in self.generate_log_spaced_embryo_cells() {
            let Some(candidate_body) = self.estimate_candidate_body(embryo_cell) else {
                continue;
            };
            let candidate_body = self.coalesce_generated_body(candidate_body);
            if candidate_body.is_trivial_mass() {
                continue;
            }
            Body::insert(&mut self.bodies, candidate_body);
        }

        self.dust_left = false;

        let log_level = *get_log_level!();
        log!(
            log_level,
            1,
            "Generated {} semi-analytic bodies between {:.2} AU and {:.2} AU",
            self.bodies.len(),
            self.planet_inner_bound,
            self.planet_outer_bound
        );
        if log_level >= 2 && !self.bodies.is_empty() {
            log!(log_level, 2, "{}", self);
        }

        self
    }

    pub fn accrete_with_mode(&mut self, generation_mode: GenerationMode) -> &mut Self {
        match generation_mode {
            GenerationMode::Aggregation => self.accrete_aggregation(),
            GenerationMode::SemiAnalyticExperimental => self.accrete_semi_analytic(),
        }
    }

    pub fn accrete(&mut self) -> &mut Self {
        self.accrete_with_mode(GenerationMode::Aggregation)
    }
}

//---------------------------  AccretionStepper  ------------------------------

/// A thin wrapper that lets the UI drive the aggregation loop one step at a
/// time.  Create it from a ready `AccretionDisk`, read the initial disk state,
/// then call `step()` repeatedly until `is_done()`.
pub struct AccretionStepper {
    disk: AccretionDisk,
}

impl AccretionStepper {
    /// Create a new stepper.  Calls `reset_generation_state()` on the disk so
    /// it is ready for the first `step()` call.
    pub fn new(mut disk: AccretionDisk) -> Self {
        disk.reset_generation_state();
        AccretionStepper { disk }
    }

    /// Returns `true` when there is no remaining dust to accrete.
    pub fn is_done(&self) -> bool {
        !self.disk.dust_left
    }

    /// Perform one aggregation step (one protoplanet injection attempt).
    pub fn step(&mut self) -> AccretionStepResult {
        self.disk.aggregation_step()
    }

    /// Borrow the underlying disk (e.g. to read `bands` or `bodies`).
    pub fn disk(&self) -> &AccretionDisk {
        &self.disk
    }

    /// Consume the stepper and return the finished disk.
    pub fn into_disk(self) -> AccretionDisk {
        self.disk
    }
}

#[cfg(test)]
mod tests {
    use super::{AccretionDisk, AccretionStepper, StepOutcome};
    use crate::{
        accretion_parameters::{get_accretion_parameters, set_accretion_parameters},
        body::Body,
        condensation::RegionThermalProfile,
        generation::GenerationMode,
        random::set_rng_seed,
        types::MassType,
    };
    use std::sync::{Mutex, MutexGuard, OnceLock};

    static TEST_ACCRETION_STATE_MUTEX: OnceLock<Mutex<()>> = OnceLock::new();

    fn lock_accretion_test_state() -> MutexGuard<'static, ()> {
        TEST_ACCRETION_STATE_MUTEX
            .get_or_init(|| Mutex::new(()))
            .lock()
            .unwrap()
    }

    fn body_at(a: f64, e: f64, mass_in_sols: f64) -> Body {
        Body {
            a,
            e,
            mass_in_sols,
            mass_type: MassType::Planet,
            ..Body::default()
        }
    }

    fn total_disk_mass_in_sols(disk: &AccretionDisk) -> f64 {
        disk.bodies.iter().map(|body| body.mass_in_sols).sum()
    }

    #[test]
    fn find_collision_detects_inner_new_body_against_outer_neighbor() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            ..AccretionDisk::default()
        };
        disk.bodies.push(body_at(1.5, 0.2, 0.01));

        let new_body = body_at(1.0, 0.0, 0.01);
        let (found_collision, closest_neighbor) = disk.find_collision(&new_body);

        assert!(found_collision);
        assert_eq!(closest_neighbor, 0);
    }

    #[test]
    fn find_collision_detects_outer_new_body_against_inner_neighbor_outward_reach() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            ..AccretionDisk::default()
        };
        disk.bodies.push(body_at(1.0, 0.2, 0.01));

        let new_body = body_at(1.5, 0.0, 0.01);
        let (found_collision, closest_neighbor) = disk.find_collision(&new_body);

        assert!(found_collision);
        assert_eq!(closest_neighbor, 0);
    }

    #[test]
    fn find_collision_detects_overpacked_ratio_below_threshold() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            ..AccretionDisk::default()
        };
        disk.bodies.push(body_at(1.0, 0.0, 0.000001));

        // 1.06 / 1.0 = 1.06 < MIN_PLANETARY_LOG_SPACING_RATIO (1.15) → overpacked via ratio check
        let new_body = body_at(1.06, 0.0, 0.00003);
        let (found_collision, closest_neighbor) = disk.find_collision(&new_body);

        assert!(found_collision);
        assert_eq!(closest_neighbor, 0);
    }

    #[test]
    fn find_collision_detects_overpacked_hill_sphere_for_massive_bodies() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            ..AccretionDisk::default()
        };
        // Bodies at 1.0 and 1.5 AU with 0.01 solar masses each pass the ratio
        // check (1.5 > 1.15) but fail the mutual Hill sphere criterion.
        disk.bodies.push(body_at(1.0, 0.0, 0.01));

        let new_body = body_at(1.5, 0.0, 0.01);
        let (found_collision, _) = disk.find_collision(&new_body);

        assert!(found_collision);
    }

    #[test]
    fn find_collision_no_collision_for_well_separated_bodies() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            ..AccretionDisk::default()
        };
        // Bodies at 1.0 and 5.0 AU with tiny masses — ratio ok (5.0), Hill ok
        disk.bodies.push(body_at(1.0, 0.0, 1e-8));

        let new_body = body_at(5.0, 0.0, 1e-8);
        let (found_collision, _) = disk.find_collision(&new_body);

        assert!(!found_collision);
    }

    #[test]
    fn coalesce_body_merges_additional_overlapping_neighbors() {
        let mut disk = AccretionDisk {
            central_mass_in_sols: 1.0,
            luminosity_in_sols: 1.0,
            thermal_profile: RegionThermalProfile::from_host_luminosity(1.0),
            ..AccretionDisk::default()
        };
        let mut empty_band = super::Band::new(0.5, 1.5);
        empty_band.dust_present = false;
        empty_band.gas_present = false;
        disk.bands.push(empty_band);
        disk.bodies.push(body_at(1.0, 0.0, 0.0001));
        disk.bodies.push(body_at(1.12, 0.0, 0.0001));

        let (merged_body, _collisions) = disk.coalesce_body(body_at(1.06, 0.0, 0.0001));

        assert!(disk.bodies.is_empty());
        assert!((merged_body.mass_in_sols - 0.0003).abs() < 1e-12);
        assert!(merged_body.a > 1.0);
        assert!(merged_body.a < 1.12);
    }

    #[test]
    fn semi_analytic_accrete_is_deterministic_for_same_seed() {
        let _guard = lock_accretion_test_state();
        let original = get_accretion_parameters();

        set_rng_seed(4242);
        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
        let mut first = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
        first.accrete_with_mode(GenerationMode::SemiAnalyticExperimental);
        let first_snapshot = first
            .bodies
            .iter()
            .map(|body| (body.a, body.e, body.mass_in_sols, body.mass_type))
            .collect::<Vec<_>>();

        set_rng_seed(4242);
        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
        let mut second = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
        second.accrete_with_mode(GenerationMode::SemiAnalyticExperimental);
        let second_snapshot = second
            .bodies
            .iter()
            .map(|body| (body.a, body.e, body.mass_in_sols, body.mass_type))
            .collect::<Vec<_>>();

        assert_eq!(first_snapshot, second_snapshot);

        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
    }

    #[test]
    fn semi_analytic_accrete_respects_planet_bounds_and_sorts_bodies() {
        let _guard = lock_accretion_test_state();
        let original = get_accretion_parameters();

        set_rng_seed(17);
        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
        let mut disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 1.0, 5.0);
        disk.accrete_with_mode(GenerationMode::SemiAnalyticExperimental);

        assert!(!disk.bodies.is_empty());
        assert!(!disk.dust_left);
        assert!(disk.bands.is_empty());
        for window in disk.bodies.windows(2) {
            assert!(window[0].a <= window[1].a);
        }
        for body in &disk.bodies {
            assert!(body.a >= 1.0);
            assert!(body.a <= 5.0);
        }

        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
    }

    #[test]
    fn semi_analytic_embryo_spacing_is_not_a_near_fixed_geometric_ladder() {
        let _guard = lock_accretion_test_state();

        set_rng_seed(23);
        let disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
        let embryo_cells = disk.generate_log_spaced_embryo_cells();

        assert!(
            embryo_cells.len() >= 5,
            "expected several embryo cells, got {}",
            embryo_cells.len()
        );

        let spacing_ratios = embryo_cells
            .windows(2)
            .map(|window| window[1].orbital_radius_au / window[0].orbital_radius_au)
            .collect::<Vec<_>>();

        let min_ratio = spacing_ratios.iter().copied().fold(f64::INFINITY, f64::min);
        let max_ratio = spacing_ratios.iter().copied().fold(0.0_f64, f64::max);

        assert!(
            max_ratio - min_ratio > 0.35,
            "expected noticeably irregular embryo spacing, got ratios {spacing_ratios:?}"
        );
    }

    #[test]
    fn aggregation_default_mode_is_sensitive_to_dust_density_for_same_seed() {
        let _guard = lock_accretion_test_state();
        let original = get_accretion_parameters();

        set_rng_seed(11);
        set_accretion_parameters(original.dust_density_coefficient * 0.5, original.percent_dust_in_cloud);
        let mut sparse_disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 30.0);
        sparse_disk.accrete();
        let sparse_mass = total_disk_mass_in_sols(&sparse_disk);
        let sparse_snapshot = sparse_disk
            .bodies
            .iter()
            .map(|body| (body.a, body.mass_in_sols, body.mass_type))
            .collect::<Vec<_>>();

        set_rng_seed(11);
        set_accretion_parameters(original.dust_density_coefficient * 2.0, original.percent_dust_in_cloud);
        let mut dense_disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 30.0);
        dense_disk.accrete();
        let dense_mass = total_disk_mass_in_sols(&dense_disk);
        let dense_snapshot = dense_disk
            .bodies
            .iter()
            .map(|body| (body.a, body.mass_in_sols, body.mass_type))
            .collect::<Vec<_>>();

        assert!(
            dense_mass > sparse_mass,
            "dense_mass={dense_mass} sparse_mass={sparse_mass}"
        );
        assert_ne!(dense_snapshot, sparse_snapshot);

        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
    }

    // ---- Stepper regression tests ----

    /// Snapshot of generated bodies for comparison.
    fn body_snapshot(disk: &AccretionDisk) -> Vec<(f64, f64, f64, MassType)> {
        disk.bodies
            .iter()
            .map(|b| (b.a, b.e, b.mass_in_sols, b.mass_type))
            .collect()
    }

    /// Run the stepper to completion and return the finished disk.
    fn run_stepper_to_completion(disk: AccretionDisk) -> AccretionDisk {
        let mut stepper = AccretionStepper::new(disk);
        while !stepper.is_done() {
            stepper.step();
        }
        stepper.into_disk()
    }

    /// Verify that the stepper path produces identical bodies to the one-shot
    /// `accrete_aggregation()` path for several seeds.
    ///
    /// This test is `#[ignore]` because it depends on the global RNG being
    /// uncontested.  Run with:
    ///
    /// ```text
    /// cargo test --lib -- --test-threads=1 --ignored
    /// ```
    #[test]
    #[ignore]
    fn stepper_vs_one_shot_equivalence() {
        let _guard = lock_accretion_test_state();
        let original = get_accretion_parameters();

        for seed in [42, 123, 999, 2024, 31415] {
            set_rng_seed(seed);
            set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
            let mut one_shot = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
            one_shot.accrete_with_mode(GenerationMode::Aggregation);
            let one_shot_snap = body_snapshot(&one_shot);

            set_rng_seed(seed);
            set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
            let stepper_disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
            let stepped = run_stepper_to_completion(stepper_disk);
            let stepped_snap = body_snapshot(&stepped);

            assert_eq!(
                one_shot_snap, stepped_snap,
                "seed {seed}: stepper produced different bodies than one-shot"
            );
        }

        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
    }

    #[test]
    fn stepper_step_result_matches_body_changes() {
        let _guard = lock_accretion_test_state();
        let original = get_accretion_parameters();

        set_rng_seed(42);
        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
        let disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
        let mut stepper = AccretionStepper::new(disk);

        while !stepper.is_done() {
            let bodies_before = stepper.disk().bodies.len();
            let result = stepper.step();
            let bodies_after = stepper.disk().bodies.len();

            match result.outcome {
                StepOutcome::BodyFormed => {
                    assert!(
                        bodies_after >= bodies_before,
                        "BodyFormed but body count did not increase (before={bodies_before}, after={bodies_after})"
                    );
                }
                StepOutcome::FailedInjection => {
                    assert_eq!(bodies_after, bodies_before, "FailedInjection but body count changed");
                }
            }

            assert!(
                result.injection_radius_au >= stepper.disk().planet_inner_bound,
                "injection_radius_au below planet_inner_bound"
            );
            assert!(
                result.injection_radius_au <= stepper.disk().planet_outer_bound,
                "injection_radius_au above planet_outer_bound"
            );
        }

        assert!(
            !stepper.disk().bands.iter().any(|b| b.dust_present()),
            "stepper done but dust bands remain"
        );

        set_accretion_parameters(original.dust_density_coefficient, original.percent_dust_in_cloud);
    }

    #[test]
    fn stepper_is_done_before_first_step_when_no_viable_disk() {
        let disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 40.0, 0.3);
        let stepper = AccretionStepper::new(disk);
        assert!(stepper.is_done());
    }

    #[test]
    fn band_accessors_match_initial_state() {
        let disk = AccretionDisk::new_with_planet_bounds(1.0, 1.0, MassType::Star, 0.3, 40.0);
        assert!(!disk.bands.is_empty());
        let band = &disk.bands[0];
        assert!(band.dust_present());
        assert!(band.gas_present());
        assert!(band.inner_edge() < band.outer_edge());
    }
}
