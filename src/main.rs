// src/main.rs

use starform_rust::{accretion_parameters::set_accretion_parameters, generate_star_system, random::set_rng_seed};
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "starform", about = "Generates a stellar system")]
struct Opts {
    // Indicates graphic output is desired
    // #[structopt(short = "g", help = "Display graphically (unimplemented)")]
    // flag_graphic: bool,

    // Indicates moons should be generated for all planets
    // #[structopt(short = "m", help = "Generate moons for each planet")]
    // flag_moons: bool,
    /// Indicates verbose output should be generated
    #[structopt(short = "v", help = "Enable verbose logging", default_value = "0")]
    loglevel: u8,

    /// Specifies the type of main star
    #[structopt(
        short = "t",
        long = "type",
        help = "Choose the spectral class (e.g. O, B, A, F, G) , spectral number (e.g. 3), luminosity id (S=supergiant, G=giant, D=white dwarf, M=main sequence) and orbit (e.g. 1 - 4). For instance G3M/3",
        default_value = ""
    )]
    star_type: String,

    /// Specifies the percent of the cloud that is dust
    #[structopt(
        short = "d",
        long = "dust",
        help = "Specify the percent of the cloud that is dust (e.g. 2 means 2%)",
        default_value = "0"
    )]
    dust_percent: f64,
}

/// Entry point for the star system simulation program.
///
/// This program initializes a star system based on command-line arguments provided by the user.
/// It leverages the `Opts` struct, which is derived from command-line options, to configure
/// the star system's properties and the program's logging level.
///
/// ## Command-line Arguments
/// - `star_type`: Specifies the type of star system to simulate, which can influence the
///   initialization parameters of the star system, such as spectral type or mass.
/// - `loglevel`: Sets the verbosity level of logging throughout the program's execution.
///
/// ## Behavior
/// - The program first parses command-line arguments to configure the logging level and determine
///   the type of star system to initialize.
/// - It then creates a `StarSystem` instance with the specified star type.
/// - Finally, it prints the details of the initialized star system to standard output.
///
/// ## Example Usage
/// Run the program with the following command line:
/// ```bash
/// cargo run -- -t "G3M/3" -v 2
/// ```
/// This will initialize a specific type of star system and set the logging level to 2.
fn main() {
    let opts = Opts::from_args();
    set_rng_seed(0); // Reset the global DETERMINISTIC_RNG with the new seed
    set_accretion_parameters(0.0, opts.dust_percent); // Update the global ACCRETION_PARAMETERS with the dust percent
    let star_system = generate_star_system(opts.loglevel, opts.star_type);

    println!("{star_system}");
}
