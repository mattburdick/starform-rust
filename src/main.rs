// src/main.rs

use crate::star_system::StarSystem;
use structopt::StructOpt;

mod accretion_disk;
mod body;
mod consts;
mod random;
mod star;
mod star_system;
mod types;

#[macro_use]
extern crate lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref LOGLEVEL: Mutex<u8> = Mutex::new(0); // Default log level
}

/// Retrieves the current global log level with error handling.
///
/// This macro accesses a globally shared mutable log level, which is stored within a Mutex to ensure thread safety.
/// It attempts to lock the Mutex and retrieve the value. If the Mutex is poisoned (i.e., a panic occurred while the
/// lock was held previously), this macro handles the error by printing a message to standard error and then
/// attempting to recover the original MutexGuard.
///
/// # Returns
/// Returns the current log level as a `u8`.
///
/// # Panics
/// The macro will panic if the Mutex is irrecoverably poisoned after an error has occurred.
///
/// # Example
/// Here's how you might use this macro in a function to adjust log verbosity based on the retrieved log level:
///
/// ```
/// if *get_log_level!() > 2 {
///     println!("Verbose logging enabled.");
/// } else {
///     println!("Verbose logging disabled.");
/// }
/// ```
///
/// # Errors
/// If an error occurs while locking the Mutex (indicating that another thread panicked while holding the lock),
/// the error is logged to stderr, and the program attempts to continue by retrieving the contained value,
/// possibly causing a second panic if the Mutex is irrecoverably damaged.
#[macro_export]
macro_rules! get_log_level {
    () => {
        $crate::LOGLEVEL.lock().unwrap_or_else(|e| {
            eprintln!("Error locking LOGLEVEL: {}", e);
            e.into_inner() // Return the MutexGuard
        })
    };
}

/// Sets the global log level to a new value in a thread-safe manner.
///
/// This macro attempts to lock a global Mutex that stores the log level and sets it to a new specified value.
/// If the Mutex is successfully locked, the log level is updated. If the Mutex cannot be locked (which may
/// happen if another thread panicked while holding the lock), an error message is printed to standard error.
///
/// # Parameters
/// - `$new_level`: The new log level to set, expressed as a `u8`.
///
/// # Usage
/// This macro is used to dynamically change the verbosity level of logging across different parts of an application
/// that share a common logging framework. It ensures that the update is thread-safe by using a Mutex to protect
/// the global log level state.
///
/// # Example
/// ```
/// // Set the global log level to 3
/// set_log_level!(3);
/// ```
///
/// # Errors
/// If the Mutex is poisoned due to a previous error (e.g., a panic in another thread while it was holding the lock),
/// this macro will output an error message to stderr indicating that the log level could not be set.
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

/// Macro to conditionally execute `println!` based on the specified log level threshold.
///
/// This macro checks if the current log level meets or exceeds the specified threshold before executing
/// the given `println!` statement. It is designed to reduce boilerplate code associated with conditional
/// logging throughout an application.
///
/// # Usage
/// - `log_level`: The current logging level of the application.
/// - `threshold`: The minimum log level required for the `println!` to execute.
/// - `fmt`: A format string as you would supply to `println!`.
/// - `args`: Comma-separated list of arguments to pass to the format string (if any).
///
/// # Example
/// ```
/// log!(current_log_level, 3,
///     "Protoplanet mass={:.2} {}, mass_density: {:.2}, non-giant density={:.2}",
///     new_mass, if collects_gas { "(gas giant)" } else { "" }, mass_density, dust_density);
/// ```
#[macro_export]
macro_rules! log {
    ($log_level:expr, $threshold:expr, $fmt:expr $(, $args:expr)*) => {
        if $log_level >= $threshold {
            println!($fmt $(, $args)*);
        }
    };
}

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
        default_value = ""
    )]
    star_type: String,
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
    set_log_level!(opts.loglevel);

    let star_system = StarSystem::new(&opts.star_type);

    println!("{}", star_system.to_string());
}
