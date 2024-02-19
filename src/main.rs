// use std::path::PathBuf;

// https://www.tenderisthebyte.com/blog/2019/05/08/parsing-cli-args-with-structopt/
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
struct Opts {
    #[structopt(short = "g", help = "Display graphically (unimplemented)")]
    flag_graphic: bool,
    #[structopt(short = "m", help = "Generate moons for each planet")]
    flag_moons: bool,
    #[structopt(short = "v", help = "Display graphically (unimplemented)")]
    flag_verbose: bool,
    #[structopt(
        short = "t",
        help = "Choose the spectral class (e.g. O, B, A, F, G) , spectral number (e.g. 3), luminosity id (S=supergiant, G=giant, D=white dwarf, M=main sequence) and orbit (e.g. 1 - 4). For instance G3M/3"
    )]
    star_type: String,
    // #[structopt(short, long, parse(from_os_str))]
    // infile: PathBuf,

    // #[structopt(parse(from_os_str))]
    // outfile: PathBuf,
}

struct Star {
    spec_class: char,
    spec_num: i32,
    lum_id: char,
}

impl Star {
    // Factory method that takes a string and returns a Result type, indicating success with `Star` or failure.
    fn from_str(input: &str) -> Result<Star, &'static str> {
        let parts: Vec<&str> = input.split('/').collect();
        if parts.len() != 2 {
            return Err("Input does not match the expected format.");
        }

        let spec_info = parts[0];
        if spec_info.len() < 3 {
            return Err("Spectral information is incomplete.");
        }

        let spec_class = spec_info.chars().next().unwrap();
        let lum_id = spec_info.chars().nth(spec_info.len() - 1).unwrap();
        let spec_num: i32 = match spec_info[1..spec_info.len() - 1].parse() {
            Ok(num) => num,
            Err(_) => return Err("Failed to parse spectral number."),
        };

        Ok(Star {
            spec_class,
            spec_num,
            lum_id,
        })
    }
}

fn main() {
    let opts = Opts::from_args();
    let star_str = &opts.star_type;
    match Star::from_str(star_str) {
        Ok(star) => println!(
            "Star created: {}{}{}",
            star.spec_class, star.spec_num, star.lum_id
        ),
        Err(e) => println!("Error creating star: {}", e),
    }
    println!("{:?}", opts);
    println!("Hello, world!");
}
