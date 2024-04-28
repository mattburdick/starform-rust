// src/accretion_disk.rs

use crate::{consts, get_log_level};
use rand::Rng;
use std::{collections::VecDeque, fmt};

use crate::body::Body;
use crate::types::MassType;

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

//---------------------------  AccretionDisk  ---------------------------------

#[derive(Debug)]
pub struct AccretionDisk {
    pub central_mass: Body,
    pub luminosity: f64,         // The luminosity of the nearby star in solar luminosities
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
    pub fn random_eccentricity() -> f64 {
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
            body.collide(injected_body);
            body.mass = self.accrete_dust(body, stell_luminosity_ratio, dust_density);

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
    fn accrete_dust(&mut self, mut body: Body, crit_mass: f64, dust_density: f64) -> f64 {
        loop {
            let old_mass = body.mass;
            body.mass = self.collect_dust(body, crit_mass, dust_density);

            if (body.mass - old_mass) <= (0.0001 * old_mass) {
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
                body.mass * consts::SUN_MASS_IN_EARTH_MASSES,
                body.a
            );
        }

        body.mass
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
    pub fn collect_dust(&mut self, protoplanet: Body, crit_mass: f64, dust_density: f64) -> f64 {
        // The injected mass has inner and outer effect limits based on size of the mass, its orbit, and eccentricity.
        // We will accumulate mass from all the dust and gas bands into "mass".
        let log_level = *get_log_level!();

        // Calculate the inner and outer effect limits of the injected mass
        // The inner and outer effect limits are the distances from the primary star at which the mass can affect the dust bands.
        // The aphelion and parahelion influences determine the distance from the orbital plane at which the mass can affect the dust bands at aphelion and parahelion.
        let (r_inner, r_outer, x_aphelion, x_parahelion) =
            Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass);

        let mut new_mass = protoplanet.mass;

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

            if log_level >= 3 {
                println!(
                    "Protoplanet mass={:.2}{}, mass_density: {:.2}, non-giant density={:.2}",
                    new_mass,
                    if collects_gas { " (gas giant)" } else { "" },
                    mass_density,
                    dust_density
                );
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

        // // Now update `dust_left` based on whether any band still has dust.
        // self.dust_left = self.bands.iter().any(|band| band.dust_present);

        new_mass
    }

    /// docs go here
    pub fn new(central_mass: Body, luminosity: f64, distance_from_primary_star_in_au: f64) -> Self {
        let planet_inner_bound = if central_mass.mass_type == MassType::Planet {
            // The inner bound for moons is the Roche limit of the planet
            let diameter_in_au = distance_from_primary_star_in_au * 2.0 / consts::KM_PER_AU;
            2.44 * diameter_in_au
        } else {
            0.3 * central_mass.mass.powf(1.0 / 3.0)
        };

        let planet_outer_bound = 50.0 * central_mass.mass.powf(1.0 / 3.0);

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
        disk_outer_bound = disk_outer_bound.min(Self::stell_dust_limit(central_mass.mass, 0.0, central_mass.mass_type));

        // Init the bands of dust / gas in the accretion disk
        let mut bands: Vec<Band> = Vec::new();
        Band::insert(&mut bands, Band::new(disk_inner_bound, disk_outer_bound));

        AccretionDisk {
            luminosity,
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
            let a: f64 = rand::thread_rng().gen_range(bound1..=bound2);
            let mut protoplanet = Body::new(a, e, consts::PROTOPLANET_MASS, MassType::Planet);
            let (eff_inner_bound, eff_outer_bound, _, _) =
                Body::gravitational_effect_limits(protoplanet.a, protoplanet.e, protoplanet.mass);

            if self.dust_available(eff_inner_bound, eff_outer_bound) {
                if log_level >= 2 {
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
                    * self.central_mass.mass.sqrt()
                    * (consts::ALPHA * -1.0 * a.powf(1.0 / consts::N)).exp();
                /*
                 *	Assume that dust is ten times more dense around planets:
                 */
                if self.central_mass.mass_type == MassType::Planet {
                    dust_density = dust_density * 10.0;
                }

                let crit_mass = Body::critical_limit(a, e, self.luminosity);
                protoplanet.mass = self.accrete_dust(protoplanet, crit_mass, dust_density);
                if (protoplanet.mass != 0.0) && (protoplanet.mass != consts::PROTOPLANET_MASS) {
                    self.coalesce_planetesimals(
                        &protoplanet,
                        crit_mass,
                        self.luminosity,
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
