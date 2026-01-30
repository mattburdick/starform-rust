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

If you want to use your library in other projects, you can add it as a dependency in the
`Cargo.toml` of the consuming project:

```toml
[dependencies]
my_project = { path = "../path_to_starform" }
```

Clone the two repos into the same folder structure. For instance:

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
