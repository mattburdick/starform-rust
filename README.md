# starform-rust

This is both a command line app and a lib.

## Building and running starform

- To build the library:

  ```sh
  cargo build --lib
  ```

- To build and run the command line app:

  ```sh
  cargo run
  ```

## Using the Library in Other Projects

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
