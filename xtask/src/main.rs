use std::process::{Command, ExitCode};

fn main() -> ExitCode {
    let mut args = std::env::args().skip(1);

    match args.next().as_deref() {
        Some("lint") => lint(),
        _ => {
            eprintln!("Usage: cargo lint (via xtask)\n\nCommands:\n  lint  Run Rust clippy and Markdown lint");
            ExitCode::from(2)
        }
    }
}

fn lint() -> ExitCode {
    let mut clippy = Command::new("cargo");
    clippy.args(["clippy", "--all-targets", "--", "-D", "warnings"]);
    if !run(clippy) {
        return ExitCode::from(1);
    }

    let mut markdownlint = Command::new("npx");
    markdownlint.args(["--yes", "markdownlint-cli2"]);
    if !run(markdownlint) {
        eprintln!(
            "\nMarkdown lint failed.\n\nIf `npx` is missing, install Node.js (which provides `npx`), then rerun `cargo lint`."
        );
        return ExitCode::from(1);
    }

    ExitCode::SUCCESS
}

fn run(mut command: Command) -> bool {
    eprintln!("+ {:?}", command);
    match command.status() {
        Ok(status) => status.success(),
        Err(err) => {
            eprintln!("Failed to run command: {err}");
            false
        }
    }
}
