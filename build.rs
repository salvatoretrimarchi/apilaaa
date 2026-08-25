// The release version does not come from Cargo.toml.
//
// CI knows what it is publishing — the tag, or one patch above the highest
// one — and used to stamp that into Cargo.toml *and* into the apilaaa entry
// of Cargo.lock, because the lockfile carries the package's own version and
// `--locked` refuses if the two disagree. Editing two generated files by
// text surgery moments before demanding they agree is one runner away from
// a build that fails for reasons that have nothing to do with the code.
//
// So nothing is stamped any more: CI sets APILAAA_VERSION and the binary
// reads it at compile time (see `VERSION` in `main.rs`). Cargo.toml and
// Cargo.lock reach the compiler exactly as they are committed, which is what
// `--locked` is supposed to check in the first place.
//
// This file exists only to tell Cargo that the variable matters, so a local
// build with a different APILAAA_VERSION is not served a stale binary.
fn main() {
    println!("cargo::rerun-if-env-changed=APILAAA_VERSION");
}
