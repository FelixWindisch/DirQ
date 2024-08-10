# DirQ: Fast Directed Q Nearness
DirQ is a Rust crate that computes directed (q, i, j)-graphs efficiently. For more information about directed Q-analysis, refer to https://arxiv.org/abs/2202.07307 or my Master's Thesis "A Novel Definition of Directed q-nearness :  Comparative Analysis and Algorithms" (soon to be published).
## Building from source
The crate has been developed for Rust and cargo version 1.68.2. Running `cargo build --release` at the project root will generate an executable in target/release.
## Running the executable
After building from source or finding the prebuilt binary in the builds directory, executing it in the command line will start the GUI. 
Alternatively, the command line interface can be used by passing the following arguments to the executable:
| Argument   |      Meaning  | Possible Values    |  Default |
|------------|:-------------:|-------------------:|----------|
| input_file | Input File in .flagser file format | File Path | None |
| output_folder | Folder to save the output .flagser file to |  File Path | None |
| q | Parameter q for (q, i, j)-nearness |  Integer | None |
| i | Parameter i for (q, i, j)-nearness |  Integer | None |
| j | Parameter j for (q, i, j)-nearness |  Integer | None |
| max_dimension | Maximal dimension of simplices in the directed flag complex |  Integer | Infinity |
| definition | Whether to use the original or the novel definition for q-nearness |  "new", "original" | "new" |
| parallelization | Method for Multithreading |  "none", "mutex", "split_and_merge", "local_cofaces", "bottom_up", "naive" | "none" |
| benchmark | Run in Benchmark mode |  "false", "parallel", "true", "full" | "false" |
