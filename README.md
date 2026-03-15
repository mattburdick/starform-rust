# starform-rust

Generate procedural star systems in Rust.

This crate is both:

- a **library** you can call from your own Rust projects, and
- a small **CLI** that prints a generated system to stdout.

If you want to *see what Starform generates* (instead of reading Rust output):

- View Starform-generated systems: [mattburdick.github.io/starform](https://mattburdick.github.io/starform/)
- Background + history: [Reintroducing Starform](https://mattburdick.github.io/reintroducing-starform/)

## What you get

Starform’s output is intentionally "classic sci-fi" in vibe: a readable report describing a star system and the bodies that form around it.
In this Rust implementation you’ll typically see:

- 1–4 star systems (single, binary, trinary, quaternary) with companion orbital radii
- Per-star accretion that builds out a set of orbiting bodies
- A text-first summary that’s easy to save, diff, or feed into other tooling

## A quick scientific note (please read)

Starform is inspired by older, classic work on planet formation and accretion models.
That means it’s great for exploration, worldbuilding, and "what if" tinkering — but it is **not a modern astrophysics simulator**.

Science has moved on since the original research paper that spawned Starform.
If you want a project that aims to reflect **more up-to-date science**, check out StarForge:
[starforge.space](https://starforge.space/index.html)

## Atmosphere modeling notes

Recent atmosphere work in `starform-rust` is being guided by a small research pass over the PDFs in
`docs/other/`, especially the papers listed below.

The main takeaway is that interesting planetary atmospheres should not be treated as a single
pressure proxy or an all-or-nothing label. A better approximation is to model atmospheres as a
retained gas mixture shaped by:

- formation location and volatile delivery
- condensation chemistry in the protoplanetary disk
- the difference between primordial `H2/He` envelopes and secondary atmospheres
- atmospheric escape and photoevaporation
- outgassing, water loss, and greenhouse/irradiation state

One concrete issue this work is addressing is false-airless planets. The older Starform-style logic
can effectively collapse a world to "airless" if it cannot retain `N2`, even when it should still be
able to hold onto heavier gases such as `CO2`, `SO2`, `H2O`, or `Ar`. The updated direction in this
crate is therefore:

- use species-by-species retention rather than a single `N2` gate
- derive pressure from retained atmospheric inventory
- expose gas percentages and hazard traits alongside compact atmosphere classes
- preserve the existing `AtmosphereClass` as a presentation layer, not the only internal model

This remains a conservative, game-friendly approximation rather than a full photochemistry solver,
but it should produce more varied results such as:

- thin `CO2` atmospheres on marginal rocky worlds
- dense `CO2` greenhouse atmospheres on hot volatile-rich planets
- `N2`-dominated temperate atmospheres with water vapor and oxygen admixtures
- hydrogen-rich mini-Neptune/sub-Neptune envelopes
- sulfurous or corrosive hot atmospheres

### Papers consulted for the atmosphere update

- `docs/other/Formation_Orbital_and_Internal_Evolution.pdf`
- `docs/other/Planetary_Dynamics_and_Habitable_Planet.pdf`
- `docs/other/Academia.edu_Bundle_-_Terrestrial_planet_formation_in_extra_so/most_similar_papers_to_this_one/Astrophysics_of_Planet_Formation.pdf`
- `docs/other/Academia.edu_Bundle_-_Terrestrial_planet_formation_in_extra_so/papers_cited_by_this_one/A_new_family_of_planets_Ocean_Planets.pdf`
- `docs/other/Academia.edu_Bundle_-_Terrestrial_planet_formation_in_extra_so/most_similar_papers_to_this_one/Chemistry_in_an_Evolving_Protoplanetary.pdf`

## Quickstart

### Build

```sh
cargo build
```

### Run (CLI)

Generate a system and print it:

```sh
cargo run
```

The CLI supports a couple useful knobs:

- Choose a star type (example from the built-in help text):

  ```sh
  cargo run -- -t "G3M/3"
  ```

- Increase verbosity:

  ```sh
  cargo run -- -v 2
  ```

- Set dust percentage for the accretion cloud:

  ```sh
  cargo run -- -d 2
  ```

Run `cargo run -- --help` to see the full set of options.

Note: the CLI currently resets the RNG seed to `0`, so runs are reproducible by default.

## Using the library

At the highest level, you can generate a system via `generate_star_system(loglevel, star_type)`.
If `star_type` is an empty string, the system is generated randomly; otherwise, you can request a specific stellar classification.

### Minimal example

```rust
use starform_rust::generate_star_system;

fn main() {
  // Empty star_type => random system
  let system = generate_star_system(0, "".to_string());
  println!("{system}");

  // Or request a specific star type
  let g_star_system = generate_star_system(0, "G3M/3".to_string());
  println!("{g_star_system}");
}
```

## Adding `starform-rust` to another project

If you want to use this library in other projects, add it as a dependency in the
`Cargo.toml` of the consuming project.

### Recommended: git dependency

This is the easiest setup for CI and for consumers who don't want a local checkout:

```toml
[dependencies]
starform-rust = { git = "https://github.com/mattburdick/starform-rust.git" }
```

### Local development: path dependency

If you're iterating on `starform-rust` itself (or on a consumer like `starform-ui`) you can use a
path dependency instead:

```toml
[dependencies]
starform-rust = { path = "../starform-rust" }
```

There are a couple common directory layouts that work well:

#### Option A: sibling repos (recommended for two-repo workflows)

```text
/Users/yourname/Documents/GitHub/
├── starform-rust/
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs
│       └── main.rs
├── starform-ui/
│   ├── Cargo.toml
│   └── src/
│       └── main.rs
```

#### Option B: vendored under a consumer repo (handy for quick local iteration)

This is also valid, but is purely a developer convenience. For example, `starform-ui` can be
configured to depend on `starform-rust = { path = "starform-rust" }` if this folder exists:

```text
starform-ui/
├── Cargo.toml
├── src/
└── starform-rust/
  ├── Cargo.toml
  └── src/
```

## References

Papers that have informed the algorithms and models in this crate. In-code
comments use short author–year citations (e.g. `// (Raymond 2008)`) that
refer back to this list.

1. **Dole, S. H.** (1969). "Formation of Planetary Systems by Aggregation: A
   Computer Simulation." *Icarus*, 13, 494–508.
   — The original band-based accretion algorithm: protoplanets are injected at
   random radii and sweep dust/gas from overlapping bands until accretion
   stalls. Constants such as the eccentricity coefficient, dust density, and
   critical mass limit originate here.

2. **Fogg, M. J.** (1985). "Extra-Solar Planetary Systems: A Microcomputer
   Simulation." *Journal of the British Interplanetary Society*, 38, 501–514.
   — Extended Dole's model with volatile inventory estimation, atmosphere
   retention, and environmental classification (used throughout `enviro.rs`).

3. **Raymond, S. N.** (2008). "Terrestrial planet formation in extra-solar
   planetary systems." *Proceedings IAU Symposium No. 249*, 233–250.
   doi:10.1017/S1743921308016645.
   — Demonstrates that terrestrial planet compositions are set by feeding-zone
   width and radial mixing during late-stage accretion, not by local
   condensation temperature alone. The feeding-zone composition tracking in
   `accretion_disk.rs` and mass-weighted composition blending during collisions
   in `body.rs` are based on this work.

4. **Fortier, A., Alibert, Y., Carron, F., Benz, W., & Dittkrist, K.-M.**
   (2013). "Planet formation models: the interplay with the planetesimal disc."
   *Astronomy & Astrophysics*, 549, A44. doi:10.1051/0004-6361/201220241.
   — Core-nucleated accretion with oligarchic growth; consulted for
   planetesimal disc dynamics and gas accretion timing.
