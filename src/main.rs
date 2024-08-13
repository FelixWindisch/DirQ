extern crate array_tool;

use dir_q::graph;

use dir_q::directed_q;
use dir_q::gui;
pub use graph::{EdgeMapGraph, EdgeListGraph, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = graph::EdgeMapGraph;



use clap::Parser;
use std::env;


//Example Usage 
//cargo run --release -- --input-file "test.flag" --i 0 --j 2 --q 2  --definition "new" --parallelization "none"
//Benchmark
//cargo run --release -- --input-file "c_elegans.flag" --i 0 --j 3 --q 3  --benchmark "true"


fn main() {
    let args_raw: Vec<String> = env::args().collect();
    if args_raw.len() < 2
    {
        let result = gui::do_gui().unwrap_or_else(|_x|{print!("UI crashed suddenly")});
        return result;
    }
    else
    {
        let args = directed_q::Args::parse();
        directed_q::run(args); 
    }    
}
