// src/accretion_disk.rs

use crate::random::get_random_number;
use crate::{body::Body, consts, get_log_level, log, types::MassType};
use std::sync::{Arc, RwLock};
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

//---------------------------  AccretionDisk  ---------------------------------

#[derive(Debug, Clone)]
pub struct AccretionDisk {
    pub central_mass_in_sols: f64,
    pub luminosity_in_sols: f64, // The luminosity of the nearby star in solar luminosities
    pub planet_inner_bound: f64, // Inner limit at which a body can exist in orbit about the central mass
    pub planet_outer_bound: f64, // Outer limit at which a body can exist in orbit about the central mass
    pub disk_inner_bound: f64,   // Inner limit of the accretion disk
    pub disk_outer_bound: f64,   // Outer limit of the accretion disk
    pub bands: Vec<Band>,        // The bands of dust and gas comprising the accretion disk
    pub bodies: Vec<Body>,       // The bodies in the accretion disk
    pub dust_left: bool,         // If true, dust remains to be accreted
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
        write!(f, "\n")?;

        for body in &self.bodies {
            let distance = body.a.round() as usize;
            let body_char = match body.mass_type {
                MassType::Moon => "o",
                MassType::Planet => "o",
                MassType::GasGiant => "0",
                MassType::Star => "*",
            };
            if distance > 0 {
                write!(f, "{}{body_char}\n", " ".repeat(distance))?;
            }
        }
        write!(f, "\n")?;

        for band in &self.bands {
            write!(f, "{}", band)?;
        }
        write!(f, "\n")?;

        Ok(())
    }
}

impl Default for AccretionDisk {
    fn default() -> Self {
        Self::new(&Body::default(), 0.0)
    }
}

impl AccretionDisk {
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
    /// ```
    /// let eccentricity = random_eccentricity();
    /// println!("Generated eccentricity: {}", eccentricity);
    /// ```
    pub fn random_eccentricity() -> f64 {
        let random_number_between_0_and_1: f64 = get_random_number(0.0001..=1.0);
        1.0 - (1.0 - random_number_between_0_and_1).powf(consts::ECCENTRICITY_COEFF)
    }

    // Calculates the outer limit of gas & dust accretion disks from a planet or star
    fn stell_dust_limit(central_mass_in_sols: f64, dist_from_primary: f64, central_mass_type: MassType) -> f64 {
        let mut outer_limit = 200.0 * central_mass_in_sols.powf(1.0 / 3.0);
        if central_mass_type == MassType::Star {
            return outer_limit;
        }

        // Otherwise the central mass is a planet
        outer_limit = outer_limit / 125.0;
        let primary_effect = dist_from_primary.powf(2.0);
        if primary_effect <= 1.0 {
            outer_limit = outer_limit * primary_effect;
        }
        outer_limit
    }

    fn dust_available(&self, inside_range: f64, outside_range: f64) -> bool {
        for band in &self.bands {
            if !band.dust_present {
                continue;
            }
            if (outside_range.min(band.outer_edge) - inside_range.max(band.inner_edge)) > 0.0 {
                return true;
            }
        }

        return false;
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
    /// ```
    /// let volume = volume(2.0, 3.0, 1.0, 1.0);
    /// println!("Volume of the cylindrical shell: {}", volume);
    /// ```
    ///
    /// This calculates the volume of a cylindrical shell with an inner radius of 2 units, an outer radius of 3 units,
    /// and a total height of 2 units.
    pub fn volume(r_inner: f64, r_outer: f64, xp: f64, xa: f64) -> f64 {
        // Total height of the cylindrical shell
        let height = xa + xp;

        // Calculating the volume
        std::f32::consts::PI as f64 * height * (r_outer.powf(2.0) - r_inner.powf(2.0))
    }

    /// Determines if a new body will collide with any existing bodies in the accretion disk.
    ///
    /// This method assesses potential collisions by comparing the orbital parameters of a new body
    /// with those of existing bodies in the system. It calculates the gravitational influence range
    /// for each body and checks if their orbits overlap, which would indicate a collision.
    ///
    /// # Parameters:
    /// - `a`: Semi-major axis of the new body in astronomical units (AU).
    /// - `e`: Eccentricity of the new body's orbit.
    ///
    /// # Returns:
    /// - `(bool, usize)`: A tuple where the first element is a boolean indicating if a collision was found,
    ///    and the second element is the index of the closest body that will collide with the new body, if any.
    ///
    /// # Process:
    /// - Iterates over all bodies in the accretion disk.
    /// - For each body, calculates the separation distance and compares it against the gravitational effect distances
    ///   (`dist1` and `dist2`), which represent the maximum reach of gravitational influence for both the new body
    ///   and the existing bodies.
    /// - Determines if the new body's orbit intersects with the orbit of any existing body by checking if the
    ///   separation is less than or equal to the gravitational influence distance of either body.
    /// - Keeps track of the closest body that the new body might collide with.
    ///
    /// # Example:
    /// ```rust
    /// let mut accretion_disk = AccretionDisk::new(...);
    /// let new_body_a = 5.0; // Semi-major axis in AU
    /// let new_body_e = 0.1; // Eccentricity
    /// let (collision_found, closest_body_index) = accretion_disk.find_collision(new_body_a, new_body_e);
    /// if collision_found {
    ///     println!("Collision is likely with body at index {}", closest_body_index);
    /// }
    /// ```
    ///
    /// # Notes:
    /// - The function returns immediately if a collision is detected, with the closest colliding body.
    /// - This approach assumes a simplified model where bodies are treated as points without physical size, focusing only on their orbital parameters.
    fn find_collision(&mut self, a: f64, e: f64) -> (bool, usize) {
        let mut closest_neighbor = 0;
        let mut found_collision = false;
        let mut closest_approach = 0.0;

        for index in 0..self.bodies.len() {
            let body = &self.bodies[index];
            let separation = body.a - a;
            /*
             *  In the following calculations, 'dist1' is the distance over
             *  which the new planet gravitationally attracts the existing
             *  planet while 'dist2' is the distance over which the existing
             *  planet affects the new one.  A collision occurs if the
             *  separation of the two planets is less than the gravitational
             *  effects distance of either.
             */
            let reduced_mass = (body.mass_in_sols / (1.0 + body.mass_in_sols)).powf(0.25);
            let dist1;
            let dist2;
            if separation > 0.0 {
                /*
                 *  The neighbor is farther from the star than our test planet:
                 */
                dist1 = (a * (1.0 + e) * (1.0 + reduced_mass)) - a;
                dist2 = body.a - (body.a * (1.0 - body.e) * (1.0 - reduced_mass));
            } else {
                /*
                 *  The new planet is farther from the star than its neighbor:
                 */
                dist1 = a - (a * (1.0 - e) * (1.0 - reduced_mass));
                dist2 = (a * (1.0 + e) * (1.0 + reduced_mass)) - a;
            };

            if (separation.abs() <= dist1.abs()) || (separation.abs() <= dist2.abs()) {
                // These two should collide.  Save the distance and check the next planet.  If it is farther away, stop looking for other collisions and return the saved one:
                if found_collision {
                    /*
                     *  Already found a planet to collide with - is this
                     *  one closer?
                     */
                    if separation.abs() < closest_approach {
                        closest_approach = separation.abs();
                        closest_neighbor = index;
                    }
                } else {
                    /*
                     *  No other planet is close enough so far, so save the
                     *  current info and continue checking:
                     */
                    closest_neighbor = index;
                    closest_approach = separation.abs();
                    found_collision = true;
                }
            }
        }

        (found_collision, closest_neighbor)
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
    /// ```
    /// let mut body = Body {
    ///     mass_in_sols: 0.05,
    ///     critical_mass_limit: 0.08,
    ///     ... // other fields
    /// };
    /// let accreted_body = accretion_disk.accrete_dust(body);
    /// println!("New mass in solar masses: {}", accreted_body.mass_in_sols);
    /// ```
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
    /// ```
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
            new_mass += volume * protoplanet.local_density(new_mass);
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
    /// ```rust
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
            planet_inner_bound = planet_inner_bound.min(0.3 * primary.mass_in_sols.powf(1.0 / 3.0));
        }

        let planet_outer_bound = 50.0 * primary.mass_in_sols.powf(1.0 / 3.0);

        /*
         *	Figure out the innermost and outermost extent of the dust/gas
         *	cloud about the object.  The dust can't be any closer to the
         *	primary than can be affected by a protoplanet with zero orbital
         *	eccentricity at the minimum distance from the primary:
         */
        let (disk_inner_bound, _, _, _) =
            Body::gravitational_effect_limits(planet_inner_bound, 0.0, consts::PROTOPLANET_MASS);
        let (_, mut disk_outer_bound, _, _) =
            Body::gravitational_effect_limits(planet_outer_bound, 0.0, consts::PROTOPLANET_MASS);

        // Depending on the type of star, the outer limit of the dust cloud may be less
        disk_outer_bound = disk_outer_bound.min(Self::stell_dust_limit(primary.mass_in_sols, 0.0, primary.mass_type));

        // Init the bands of dust / gas in the accretion disk
        let mut bands: Vec<Band> = Vec::new();
        Band::insert(&mut bands, Band::new(disk_inner_bound, disk_outer_bound));

        AccretionDisk {
            luminosity_in_sols,
            central_mass_in_sols: primary.mass_in_sols,
            planet_inner_bound,
            planet_outer_bound,
            disk_inner_bound,
            disk_outer_bound,
            bands,
            bodies: Vec::new(),
            dust_left: true,
        }
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
    /// ```rust
    /// let mut accretion_disk = AccretionDisk::new(...);
    /// accretion_disk.accrete();
    /// println!("Accretion process completed: {}", accretion_disk);
    /// ```
    ///
    /// # Panics
    /// Panics if no dust band with dust is found, indicating an error in managing the state of the disk.
    /// Also, it will panic if the calculated bounds for the planet's orbit are out of expected range,
    /// which could indicate incorrect initialization or corruption of disk state.
    pub fn accrete(&mut self) -> &mut Self {
        while self.dust_left {
            let e = Self::random_eccentricity();
            let log_level = *get_log_level!();

            // Find the first band with dust
            let band = self
                .bands
                .iter()
                .find(|band| band.dust_present)
                .expect("AccretionDisk.accrete: unable to find a band with dust");

            /*
             *	Choose a location for the proto-mass that is somewhere
             *	within gravitational effect range of the first band.  As
             *	this is done for each band, the innermost band with dust
             *	still remaining will move further and further from the primary
             *	until all dust in the system has been accreted.
             */

            let bound1 = band
                .inner_edge
                .max(self.planet_inner_bound)
                .min(self.planet_outer_bound);
            let bound2 = band
                .outer_edge
                .min(self.planet_outer_bound)
                .max(self.planet_inner_bound);
            if self.planet_inner_bound > self.planet_outer_bound {
                panic!("ERROR: planet inner bound is greater than outer bound\n");
            }
            let a: f64 = get_random_number(bound1..=bound2);

            // Create a protoplanet annotated with the dust density and critical mass limit at this distance from the primary
            let mut protoplanet = Body::new(
                a,
                e,
                consts::PROTOPLANET_MASS,
                MassType::Planet,
                0.0, // The radius of the protoplanet is not used
                self.central_mass_in_sols,
                self.luminosity_in_sols,
                Some(Arc::new(RwLock::new(AccretionDisk::default()))),
            );
            let (eff_inner_bound, eff_outer_bound, _, _) =
                Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass_in_sols);

            if !self.dust_available(eff_inner_bound, eff_outer_bound) {
                // This shouldn't happen as we've already checked for dust within gravitational effect range
                panic!("ERROR: found no dust in range of protoplanet");
            }
            log!(
                log_level,
                2,
                "    Injecting proto-{} at {:.2} AU.",
                protoplanet.mass_type,
                protoplanet.a
            );

            protoplanet = self.accrete_dust(protoplanet);
            if protoplanet.is_trivial_mass() {
                log!(
                    log_level,
                    2,
                    "    Proto-{} has trivial mass ({:.3} Earth masses).",
                    protoplanet.mass_type,
                    protoplanet.mass_in_earth_masses()
                );
                continue;
            }

            // Check for collisions with other planets
            let (found_collision, closest_neighbor) = self.find_collision(protoplanet.a, protoplanet.e);
            if found_collision {
                // We must temporarily take ownership of the body to avoid borrowing issues.
                let mut body = std::mem::replace(&mut self.bodies[closest_neighbor], Body::default());
                body.collide(&protoplanet);

                // Since it has grown in size, we need to check for more matter to accrete
                // body.mass_in_sols = self.accrete_dust(body.clone());
                body = self.accrete_dust(body.clone());

                // Now replace the original body back with the updated body.
                self.bodies[closest_neighbor] = body;
            } else {
                // The new planet won't collide with any other planet or star, so add it to the system
                if protoplanet.mass_in_sols >= protoplanet.critical_mass_limit {
                    protoplanet.mass_type = MassType::GasGiant;
                }

                protoplanet.initialize(self.luminosity_in_sols);
                Body::insert(&mut self.bodies, protoplanet);
            }

            // Create a "diagram" of the system
            log!(log_level, 1, "{}", self);
        }
        self
    }
}
