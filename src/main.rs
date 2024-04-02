pub mod consts;
pub mod random;
mod star;
pub mod star_system;
use crate::star_system::StarSystem;

// https://www.tenderisthebyte.com/blog/2019/05/08/parsing-cli-args-with-structopt/
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
    flag_verbose: u8,

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

    // Use from_str if the -t option was used on the command line. Otherwise randomly generate
    // the system
    let star_system: StarSystem;
    if &opts.star_type != "random" {
        match StarSystem::from_str(&opts.star_type) {
            Ok(value) => star_system = value,
            Err(err) => panic!("Error: {}", err),
        }
    } else {
        star_system = StarSystem::random();
    }

    // if opts.flag_verbose >= 1 {
    //     println!("Creating system with {} stars.", star_system.star_count);
    // }

    println!("{}", star_system.to_string());

    // println!("{:?}", opts);
}
