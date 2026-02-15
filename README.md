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
