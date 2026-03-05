#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: cargo (Rust) is required" >&2
  exit 2
fi

echo "+ cargo fmt --all -- --check"
cargo fmt --all -- --check

if [[ -f "Cargo.lock" ]]; then
  echo "+ cargo clippy --all-targets --all-features --locked -- -D warnings"
  cargo clippy --all-targets --all-features --locked -- -D warnings
else
  echo "+ cargo clippy --all-targets --all-features -- -D warnings"
  cargo clippy --all-targets --all-features -- -D warnings
fi

echo "ok: CI-like checks passed"
