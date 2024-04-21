pub mod accretion_disk;
pub mod consts;
pub mod random;
mod star;
pub mod star_system;
use crate::star_system::StarSystem;

#[macro_use]
extern crate lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref LOGLEVEL: Mutex<u8> = Mutex::new(0); // Default log level
}

/// Macro to get the log level safely with error handling
#[macro_export]
macro_rules! get_log_level {
    () => {
        $crate::LOGLEVEL.lock().unwrap_or_else(|e| {
            eprintln!("Error locking LOGLEVEL: {}", e);
            e.into_inner() // Return the MutexGuard
        })
    };
}
/// Macro to set the log level safely
#[macro_export]
macro_rules! set_log_level {
    ($new_level:expr) => {
        if let Ok(mut log_level) = LOGLEVEL.lock() {
            *log_level = $new_level;
        } else {
            eprintln!("Failed to lock LOGLEVEL for writing.");
        }
    };
}

use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "starform", about = "Generates a stellar system")]
struct Opts {
    /// Indicates graphic output is desired
    // #[structopt(short = "g", help = "Display graphically (unimplemented)")]
    // flag_graphic: bool,

    /// Indicates moons should be generated for all planets
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
        default_value = "random"
    )]
    star_type: String,
}

fn main() {
    let opts = Opts::from_args();
    set_log_level!(opts.loglevel);

    // If "-t" was used to specify the star type (e.g. G3M/1), generate the requested star. Otherwise randomly generate the system
    let star_system: StarSystem;
    if &opts.star_type != "random" {
        match StarSystem::from_str(&opts.star_type) {
            Ok(value) => star_system = value,
            Err(err) => panic!("Error: {}", err),
        }
    } else {
        star_system = StarSystem::random();
    }

    println!("{}", star_system.to_string());
}
