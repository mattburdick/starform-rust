use crate::{consts, get_log_level};
use rand::Rng;
use std::{collections::VecDeque, fmt};

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum MassType {
    Star,
    Planet,
    GasGiant,
}

//---------------------------  CentralMass  ------------------------------------
#[derive(Debug)]
pub struct CentralMass {
    pub mass_in_sols: f64,
    pub mass_type: MassType,
    pub luminosity_in_sols: f64,
    pub radius_in_au: f64,
}

impl CentralMass {
    pub fn new(mass_type: MassType, mass: f64, luminosity: f64, radius: f64) -> Self {
        CentralMass {
            mass_in_sols: mass,
            mass_type: mass_type,
            luminosity_in_sols: luminosity,
            radius_in_au: radius,
        }
    }
}

//---------------------------  Band  ------------------------------------------
// This represents a band of dust/gas in the accretion disk
#[derive(Debug)]
pub struct Band {
    dust_present: bool,
    gas_present: bool,
    inner_edge: f64,
    outer_edge: f64,
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

//---------------------------  Body  ------------------------------------------
#[derive(Debug)]
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
    pub fn new(a: f64, e: f64, mass: f64, mass_type: MassType) -> Self {
        Body { a, e, mass, mass_type }
    }

    fn insert(bodies: &mut Vec<Body>, new_body: Body) {
        // Find the first position where 'a' is greater
        let index = bodies.iter().position(|x| new_body.a < x.a).unwrap_or(bodies.len());

        bodies.insert(index, new_body);
    }

    fn gravitational_effect_limits(a: f64, e: f64, mass: f64) -> (f64, f64, f64, f64) {
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

    fn collide_planets(&mut self, other: &Body) {
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

        /*
         *  If the protoplanet had the misfortune to collide with a
         *  star, update the corresponding star node for the node
         *  planet node:
         */
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

//---------------------------  AccretionDisk  ---------------------------------

#[derive(Debug)]
pub struct AccretionDisk {
    pub central_mass: CentralMass,
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
        write!(f, "{:?}\n", self)
    }
}

impl AccretionDisk {
    fn random_eccentricity() -> f64 {
        let random_num: f64 = rand::thread_rng().gen_range(0.0001..=1.0);
        1.0 - random_num.powf(consts::ECCENTRICITY_COEFF)
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
            if band.outer_edge < inside_range && band.dust_present {
                return true;
            } else if band.inner_edge < outside_range && band.dust_present {
                return true;
            }
        }

        return false;
    }

    /*--------------------------------------------------------------------------*/
    /*   Orbital radius is in AU, eccentricity is unitless, and the stellar     */
    /*  luminosity ratio is with respect to the sun.  The value returned is the */
    /*  mass at which the planet begins to accrete gas as well as dust, and is  */
    /*  in units of solar masses.                                               */
    /*--------------------------------------------------------------------------*/
    fn critical_limit(orb_radius: f64, eccentricity: f64, stell_luminosity_ratio: f64) -> f64 {
        let perihelion_dist = orb_radius - orb_radius * eccentricity;
        consts::B * (perihelion_dist * stell_luminosity_ratio.sqrt()).powf(-0.75)
    }

    /*--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*/
    /*--------------------------------------------------------------------------*/
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
            let reduced_mass = (body.mass / (1.0 + body.mass)).powf(0.25);
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

    /*--------------------------------------------------------------------------*/
    /*  coalesce_planetesimals checks if the protoplanet described by a, e,     */
    /*  mass, etc, crosses the orbits of any planets already generated.  If so, */
    /*  the masses of the two planets are added (therefore assuming a perfectly */
    /*  inelastic collision) and an orbit for the resulting new planet          */
    /*  computed.  The new mass is then allowed to accrete more dust.           */
    /*  If no collision with another planet occurrs, a new planet is created    */
    /*  its statistics filled in with those of the protoplanet.                 */
    /*--------------------------------------------------------------------------*/
    fn coalesce_planetesimals(
        &mut self,
        injected_body: &Body,
        crit_mass: f64,
        stell_luminosity_ratio: f64,
        mut mass_type: MassType,
        dust_density: f64,
    ) {
        let log_level = *get_log_level!();
        if injected_body.mass <= consts::TRIVIAL_MASS {
            if log_level >= 1 {
                println!(
                    "  Trivial mass ({:.3} Earth masses) - not adding it.\n",
                    injected_body.mass * consts::SUN_MASS_IN_EARTH_MASSES
                );
            }
            return;
        }

        let (found_collision, closest_neighbor) = self.find_collision(injected_body.a, injected_body.e);
        if found_collision {
            // We must temporarily take ownership of the body to avoid borrowing issues.
            let mut body = std::mem::replace(&mut self.bodies[closest_neighbor], Body::default());
            body.collide_planets(injected_body);
            body.mass = self.accrete_dust(body.a, body.e, body.mass, stell_luminosity_ratio, dust_density);

            // Now, replace the original body back with the updated body.
            self.bodies[closest_neighbor] = body;
        } else {
            /*
             *  The new planet won't collide with any other planet or star,
             *  so allocate space for it and insert it into the system's
             *  linked list:
             */
            if injected_body.mass >= crit_mass {
                mass_type = MassType::GasGiant;
            }
            let body = Body::new(injected_body.a, injected_body.e, injected_body.mass, mass_type);
            if log_level >= 3 {
                println!("      Creating a new planet.\n");
            }

            Body::insert(&mut self.bodies, body);
        }
    }

    /*--------------------------------------------------------------------------*/
    /*  Given a mass at a particular orbit, this function repeatedly calls      */
    /*  'collect_dust' to sweep up any dust and gas it can.  Each successive    */
    /*  call to 'collect_dust' is done with the original mass plus additional   */
    /*  mass from sweeping up dust previously.  The process stops when the mass */
    /*  accumulation slows.                                                     */
    /*--------------------------------------------------------------------------*/
    fn accrete_dust(&mut self, a: f64, e: f64, mass: f64, crit_mass: f64, dust_density: f64) -> f64 {
        let mut new_mass = mass;
        loop {
            let old_mass = new_mass;
            new_mass = self.collect_dust(a, e, new_mass, crit_mass, dust_density);

            if (new_mass - old_mass) <= (0.0001 * old_mass) {
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

        if *get_log_level!() >= 1 {
            println!(
                "  Built a planet of {:0.3} Earth masses at {:0.3} AU\n",
                new_mass * consts::SUN_MASS_IN_EARTH_MASSES,
                a
            );
        }

        new_mass
    }

    /*--------------------------------------------------------------------------*/
    /*  This routine compares the location of a test mass with the location     */
    /*  of any dust and gas bands remaining.  Any dust band that lies within    */
    /*  the object's range of gravitational effect (from 'r_inner' to 'r_outer')*/
    /*  is swept up by the object.  Additionally, if the object's mass is       */
    /*  greater than the critical mass, gas is also collected by the object.    */
    /*                                                                          */
    /*  The new mass for the object is returned.                                */
    /*                                                                          */
    /*  VARIABLES PASSED IN:                                                    */
    /*    mass         Mass of the accreting object (Solar masses)              */
    /*    a            Distance from the primary star (AUs)                     */
    /*    e            Eccentricity of the object's orbit                       */
    /*    crit_mass    Mass at which, for this orbit and star, a normal planet  */
    /*                 begins to sweep up gas as well as dust and become a      */
    /*                 gas giant.                                               */
    /*    dust_head    Pointer to the head of the dust band list                */
    /* */
    /*  LOCAL VARIABLES:                                                        */
    /*    r_inner      Innermost gravitational effect limit of the object       */
    /*    r_outer      Outermost gravitational effect limit of the object       */
    /*                                                                          */
    /*--------------------------------------------------------------------------*/
    pub fn collect_dust(&mut self, a: f64, e: f64, mass: f64, crit_mass: f64, dust_density: f64) -> f64 {
        // The injected mass has inner and outer effect limits based on size of the mass, its orbit, and eccentricity.
        // We will accumulate mass from all the dust and gas bands into "mass".
        //let mut mass = injected_mass;

        // Calculate the inner and outer effect limits of the injected mass
        // The inner and outer effect limits are the distances from the primary star at which the mass can affect the dust bands.
        // The aphelion and parahelion influences determine the distance from the orbital plane at which the mass can affect the dust bands at aphelion and parahelion.
        let (r_inner, r_outer, x_aphelion, x_parahelion) = Body::gravitational_effect_limits(a, e, mass);

        let mut new_mass = mass;

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

            let collects_gas = new_mass >= crit_mass;
            if (!collects_gas && !band.dust_present) || (collects_gas && !band.gas_present && !band.dust_present) {
                // The mass is too small to pick up gas and the band has no dust.
                continue;
            }

            let mass_density = if !collects_gas {
                dust_density
            } else {
                // It's a gas giant, so the density of the dust is reduced.
                // TODO: I'm not sure what this formula is doing.
                consts::K * dust_density / (1.0 + (crit_mass / new_mass).sqrt() * (consts::K - 1.0))
            };
            print!(
                "Protoplanet mass={:.2}{}, mass_density: {:.2}, non-giant density={:.2}\n",
                new_mass,
                if collects_gas { " (gas giant)" } else { "" },
                mass_density,
                dust_density
            );

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
            band.gas_present = if collects_gas { false } else { band.gas_present };
            let swept_bandwidth = band.outer_edge - band.inner_edge;

            let volume = 2.0 * consts::PI * swept_bandwidth * (x_aphelion + x_parahelion);

            new_mass += volume * mass_density;
        }

        // Apply all changes here
        for band in new_bands {
            Band::insert(&mut self.bands, band);
        }

        Band::merge_neighbors(&mut self.bands);

        // /*
        //  *  Find the effective mass and its range of effect ('r_inner'
        //  *  through 'r_outer'):
        //  */
        // let mut accumulated_mass = mass;
        // let reduced_mass = (mass / (1.0 + mass)).powf(0.25);
        // let r_inner = Self::inner_effect_limit(a, e, reduced_mass);
        // let r_outer = Self::outer_effect_limit(a, e, reduced_mass);

        // if r_inner < 0.0 {
        //     panic!("Inner effect limit calculation error!");
        // }

        // /*
        //  *  Visit each dust band and add up any dust collected from each.
        //  *  Start with the original mass of the object:
        //  */
        // for band in self.bands.iter_mut() {
        //     /*
        //      *  If there is no gas in this band OR if the band lies outside
        //      *  the range of effect completely OR if no dust is present and
        //      *  the mass is too small to pick up gas, go to the next band:
        //      */
        //     if !band.gas_present {
        //         continue;
        //     }
        //     if band.outer_edge <= r_inner || band.inner_edge >= r_outer {
        //         continue;
        //     }
        //     if mass < crit_mass && !band.dust_present {
        //         continue;
        //     }

        //     /*
        //      *  Dust or gas exists and lies within range - sweep some up:
        //      */
        //     let mass_density = if mass < crit_mass {
        //         dust_density
        //     } else {
        //         consts::K * dust_density / (1.0 + (crit_mass / mass).sqrt() * (consts::K - 1.0))
        //     };

        //     let bandwidth = r_outer - r_inner;
        //     let temp1 = (r_outer - band.outer_edge).max(0.0);
        //     let temp2 = (band.inner_edge - r_inner).max(0.0);
        //     let width = bandwidth - temp1 - temp2;
        //     let volume = 4.0
        //         * consts::PI
        //         * a.powi(2)
        //         * reduced_mass
        //         * (1.0 - e * (temp1 - temp2) / bandwidth)
        //         * width;
        //     accumulated_mass += volume * mass_density;
        // }

        // /*
        //  *	Now re-visit each band removing dust and gas if necessary.  It
        //  *	may also be necessary to reduce the size of some dust bands
        //  *	and create new gas bands (if the dust is removed from a dust
        //  *	band, it becomes a gas band).
        //  */
        // let mut i = 0;
        // while i < self.bands.len() {
        //     let band = &mut self.bands[i];
        //     if !band.gas_present {
        //         continue;
        //     }
        //     if band.outer_edge <= r_inner || band.inner_edge >= r_outer {
        //         continue;
        //     }
        //     if mass < crit_mass && !band.dust_present {
        //         continue;
        //     }

        //     let temp1 = (r_outer - band.outer_edge).max(0.0);
        //     let temp2 = (band.inner_edge - r_inner).max(0.0);

        //     /*
        //      *  Some dust has been swept up, so update this band:
        //      */
        //     /*
        //      *  Case 1: the area of effect lies entirely within the dust band:
        //      *  Result: divide the original dust band into two smaller ones:
        //      */
        //     if temp1 == 0.0 && temp2 == 0.0 {
        //         let mut newband = Band::new(r_outer, band.outer_edge);
        //         newband.dust_present = band.dust_present;
        //         newband.gas_present = band.gas_present;
        //         self.bands.insert(i + 1, newband);
        //         band.outer_edge = r_inner;
        //         if (self.loglevel > 3) {
        //             println!(
        //                 "      Creating a new dust band 1 between ({:.2} AU - {:.2} AU).",
        //                 newband.inner_edge, newband.outer_edge
        //             );
        //         }
        //         if mass < crit_mass {
        //             /*
        //              *  The mass isn't a gas giant, so it'll sweep away all the
        //              *  dust in it's range, but leave the gas.  Therefore, we
        //              *  need to create a new gas band here:
        //              */
        //             let mut gasband = Band::new(r_inner, r_outer);
        //             gasband.dust_present = false;
        //             gasband.gas_present = true;
        //             self.bands.insert(i + 1, gasband);
        //             band.outer_edge = r_inner;

        //             if (self.loglevel > 3) {
        //                 println!(
        //                     "      Creating a new gas band 2 ({:.2} AU - {:.2} AU).",
        //                     gasband.inner_edge, gasband.outer_edge
        //                 );
        //             }
        //         }
        //     }

        //             /*
        //  *  Case 2: the area of effect encompasses the dust band entirely:
        //  *  Result: Remove the band if both dust and gas can be swept up,
        //  *          otherwise, just remove the dust from it:
        //  */
        //     else if temp1 > 0.0 && temp2 > 0.0 {
        //         if mass >= crit_mass {
        //             self.bands.remove(i);
        //             if self.loglevel > 3 {
        //                 println!(
        //                     "      Freeing a gas band 3 ({:.2} AU - {:.2} AU).",
        //                     band.inner_edge, band.outer_edge
        //                 );
        //             }
        //         } else {
        //             band.dust_present = false;
        //             if self.loglevel > 3 {
        //                 println!(
        //                     "      Removing dust from a dust/gas band 4 ({:.2} AU - {:.2} AU).",
        //                     band.inner_edge, band.outer_edge
        //                 );
        //             }
        //         }
        //     }

        //             /*
        //  *  Case 3: the area of effect and the dust band overlap with
        //  *          the dust band slightly further from the primary star.
        //  *  Result: Remove the inner part of the band if both dust and gas
        //  *          can be swept up, otherwise, just remove the dust from it:
        //  */
        //  else if temp2 > 0.0 {

        //     if mass >= crit_mass {
        //         band.inner_edge = r_outer;
        //         if self.loglevel >= 3 {
        //             println!("      Reducing a gas band 5 ({:.2} - {:.2}).",
        //                      band.inner_edge, band.outer_edge);
        //         }
        //     }
        //     else {
        //         /*
        //          *  If there is a dust band prior to this one in the list
        //          *  and it only has gas in it already and touches the inner
        //          *  edge of the current band, we don't need to create a
        //          *  new band - just add the current one onto the prior one.
        //          */
        //                         prev_band = prior_dust_band(dust_head, band);
        //                         if prev_band != NULL {
        //                             if !prev_band.dust_present
        //                                 && prev_band.outer_edge == band.inner_edge {
        //                                 prev_band.outer_edge = r_outer;
        //                                 band.inner_edge = r_outer;
        //                                 if self.loglevel >= 3 {
        //                                 	println!("      Increasing a gas band 6 ({:.2} - {:.2}).",
        //         	                        	prev_band.inner_edge, prev_band.outer_edge);
        //                                 	println!("      Reducing a dust band 7 ({:.2} - {:.2}).",
        //         	                        	band.inner_edge, band.outer_edge);
        //                                 }
        //                                 continue;
        //                             }
        //                         }
        //                         let mut gasband = Band::new(band.inner_edge, r_outer);
        //                         gasband.dust_present = false;
        //                         gasband.gas_present = true;
        //                         band.inner_edge = r_outer;
        //                         if self.loglevel >= 3 {
        //                             println!("      Reducing a dust band 8 ({:.2} - {:.2}).",
        //                                 band.inner_edge, band.outer_edge);
        //                             println!("      Creating a new gas band 9 ({:.2} - {:.2}).",
        //                                 gasband.inner_edge, gasband.outer_edge);
        //                         }

        //     }

        //             /*
        //  *  Case 4: the area of effect and the dust band overlap with
        //  *          the dust band slightly closer to the primary star.
        //  *  Result: Remove the outer part of the band if both dust and gas
        //  *          can be swept up, otherwise, just remove the dust from it:
        //  */
        //      else {
        //         if temp1 > 0.0 {
        //             if mass >= crit_mass {
        //                 band.outer_edge = r_inner;
        //             } else {
        //                                 /*
        //          *  As above, if there is a dust band after this one
        //          *  and it only has gas in it already and touches the outer
        //          *  edge of the current band, we don't need to create a
        //          *  new band - just add the current one onto the next one.
        //          */
        //             }

        //  }

        // }

        // // Now update `dust_left` based on whether any band still has dust.
        // self.dust_left = self.bands.iter().any(|band| band.dust_present);

        new_mass
    }

    /// docs go here
    pub fn new(central_mass: CentralMass, distance_from_primary_star_in_au: f64) -> Self {
        let planet_inner_bound = if central_mass.mass_type == MassType::Planet {
            // The inner bound for moons is the Roche limit of the planet
            let diameter_in_au = distance_from_primary_star_in_au * 2.0 / consts::KM_PER_AU;
            2.44 * diameter_in_au
        } else {
            0.3 * central_mass.mass_in_sols.powf(1.0 / 3.0)
        };

        let planet_outer_bound = 50.0 * central_mass.mass_in_sols.powf(1.0 / 3.0);

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
        disk_outer_bound = disk_outer_bound.min(Self::stell_dust_limit(
            central_mass.mass_in_sols,
            0.0,
            central_mass.mass_type,
        ));

        // Init the bands of dust / gas in the accretion disk
        let mut bands: Vec<Band> = Vec::new();
        Band::insert(&mut bands, Band::new(disk_inner_bound, disk_outer_bound));

        AccretionDisk {
            central_mass,
            planet_inner_bound,
            planet_outer_bound,
            disk_inner_bound,
            disk_outer_bound,
            bands,
            bodies: Vec::new(),
            dust_left: true,
        }
    }

    pub fn accrete(&mut self) -> &mut Self {
        // return self;
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

            let bound1 = band.inner_edge.max(self.planet_inner_bound).min(self.planet_outer_bound);
            let bound2 = band.outer_edge.min(self.planet_outer_bound).max(self.planet_inner_bound);
            if self.planet_inner_bound > self.planet_outer_bound {
                panic!("ERROR: planet inner bound is greater than outer bound\n");
            }
            let a: f64 = rand::thread_rng().gen_range(bound1..=bound2);
            let mut protoplanet = Body::new(a, e, consts::PROTOPLANET_MASS, MassType::Planet);
            let (eff_inner_bound, eff_outer_bound, _, _) =
                Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass);

            if self.dust_available(eff_inner_bound, eff_outer_bound) {
                if log_level >= 1 {
                    println!(
                        "  Injecting proto-{} ({:4.2} AU)",
                        if self.central_mass.mass_type == MassType::Star {
                            "planet"
                        } else {
                            "moon"
                        },
                        a
                    );
                }

                let mut dust_density = consts::DUST_DENSITY_COEFF
                    * self.central_mass.mass_in_sols.sqrt()
                    * (consts::ALPHA * -1.0 * a.powf(1.0 / consts::N)).exp();
                /*
                 *	Assume that dust is ten times more dense around planets:
                 */
                if self.central_mass.mass_type == MassType::Planet {
                    dust_density = dust_density * 10.0;
                }

                let crit_mass = Self::critical_limit(a, e, self.central_mass.luminosity_in_sols);
                protoplanet.mass =
                    self.accrete_dust(protoplanet.a, protoplanet.e, protoplanet.mass, crit_mass, dust_density);
                if (protoplanet.mass != 0.0) && (protoplanet.mass != consts::PROTOPLANET_MASS) {
                    self.coalesce_planetesimals(
                        &protoplanet,
                        crit_mass,
                        self.central_mass.luminosity_in_sols,
                        self.central_mass.mass_type,
                        dust_density,
                    );
                } else if log_level >= 2 {
                    // TODO: what does this mean?
                    println!("    Neighbor too near ({:.2} AU).", a);
                }
            } else if log_level >= 2 {
                println!("    Not enough dust at {a} AU.\n");
            }
        }
        self
    }
}
