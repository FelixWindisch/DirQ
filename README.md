# DirQ: Fast Directed Q Nearness
DirQ is a Rust crate that computes directed (q, i, j)-graphs efficiently. For more information about directed Q-analysis, refer to https://arxiv.org/abs/2202.07307 or my Master's Thesis "A Novel Definition of Directed q-nearness :  Comparative Analysis and Algorithms".
## Building from source
The crate has been developed for Rust and cargo version 1.68.2. Running `cargo build --release` at the project root will generate an executable in target/release.
## Running the executable
After building from source or finding the prebuilt binary in the builds directory, executing it in the command line will start the GUI. 
Alternatively, the command line interface can be used by passing --gui "false" to the binary. 
