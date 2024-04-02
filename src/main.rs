extern crate array_tool;
use dir_q_lib::complex;
use dir_q_lib::graph;
use dir_q_lib::graph_io;
use dir_q_lib::directed_q;
pub use graph::{EdgeMapGraph, EdgeListGraph, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = graph::EdgeMapGraph;
use std::time::Instant;
use std::collections::HashMap;
use core::iter::zip;




//Example Usage 
//cargo run --release -- --input-file "test.flag" --i 0 --j 2 --q 2  --definition "new" --parallelization "naive"
//Benchmark
//cargo run --release -- --input-file "c_elegans.flag" --i 0 --j 3 --q 3 --parallelization "naive" --benchmark "true"


use clap::Parser;
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
    // Configuration
    /// flag input file location
    #[clap(long, default_value="test.flag")]
    input_file: String,

//    /// TODO: q,i,j as a triplet
//    /// direction
//    #[clap(short, long, value_names=&["q", "i", "j"], default_value=vec![3,0,999])]
//    direction: Vec<Node>,

    /// direction_from
    #[clap(short, long, default_value="0")]
    i: usize,

    /// direction_to
    #[clap(short, long, default_value="2")]
    j: usize,

    /// q
    #[clap(short, long, default_value="2")]
    q: usize,

    /// maximal simplex dimension
    #[clap(short, long, default_value="10000000")]
    max_dimension: usize,

    /// definition: old (Henries) or new (Florians/Felix')
    #[clap(short, long, possible_values=&["old", "new"], default_value="old")]
    definition: String,

    // Sampler configuration: Technical
    /// short human-readable label for reference
    #[clap(short, long, default_value="0")]
    label: String,

    /// parallelization method
    #[clap(short, long, possible_values=&["none", "mutex", "split_and_merge", "local_cofaces", "bottom_up", "naive"], default_value="bottom_up")]
    parallelization: String,

    #[clap(short, long, possible_values=&["false", "true", "full"], default_value="false")]
    benchmark: String
}



fn main() {
    let args = Args::parse();
    let i = args.i;
    let j = args.j;
    let q = args.q;
    let max_dimension = args.max_dimension;
    assert!(max_dimension > q);
    let par = args.parallelization;
    let definition = args.definition;
    let bm = args.benchmark;
    let input_file  = args.input_file;
    let g: EdgeMapGraph = if input_file == "random.flag"
    {
        EdgeMapGraph::gen_seo_er(278, 0.028491286 *2.0, &mut rand::thread_rng())
    }
    else
    {
         graph_io::read_flag_file(&input_file)
    };
    dbg!(g.edges().len());
    if bm == "true" || bm == "full"
    {
        benchmark(g, q, i, j, bm);
        return;
    } 
    if par == "bottom_up"
    {
        if definition == "old"
        {
            let result = directed_q::old::bottom_up::get_q_digraph(&g, q, i, j);
            dbg!(result.len());
        }
        else
        {
            let result = directed_q::new::bottom_up::get_q_digraph(&g, q, i, j);
            dbg!(result.len());
        }
            
        
        return;
    }

    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);
    dbg!("Flag complex computed");
    assert!(flag_complex.len() > 1);
    if flag_complex[q+1].len() == 0
    {
        panic!("There are no q+1 simplices in the graph");
    }
    

    let q:EdgeListGraph = match  definition.as_str() 
    {
        "new" => 
        {
            match par.as_str()
            {
                "naive" => {directed_q::new::naive::get_q_digraph( g ,q, i, j, &flag_complex)}
                "mutex" => {directed_q::new::mutex::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "split_and_merge" => {directed_q::new::split_and_merge::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "local_cofaces" => {directed_q::new::local_cofaces::get_q_digraph(&g, q, i, j, &flag_complex, &simplex_map)}
                _ => {directed_q::new::single_thread::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
            }
        }
        "old" =>
        {
            match par.as_str()
            {
                "naive" => {directed_q::old::naive::get_q_digraph(g, q, i, j, &flag_complex)}
                "mutex" => {directed_q::old::mutex::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "split_and_merge" => {directed_q::old::split_and_merge::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "local_cofaces" => {directed_q::old::local_cofaces::get_q_digraph( &g, q, i, j, &flag_complex, &simplex_map)}
                _ => {directed_q::old::single_thread::get_q_digraph( q, i, j, &flag_complex,  &simplex_map)}
            }
        }
        _ => {directed_q::new::single_thread::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
    };
    let name = input_file.split('.').next().unwrap();
    graph_io::save_flag_file(&(name.to_owned() + &("_out.flag")), &q);

    //dbg!(q_graph_henri.edges().len());
    //dbg!(q_graph.edges().len());
    
    
}


fn benchmark(g: EdgeMapGraph, q:usize, i:usize, j:usize, mode: String)
{
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);

    let new_functions: Vec<for<'a, 'b> fn(usize, usize, usize, &'a Vec<Vec<Vec<u64>>>, &'b HashMap<Vec<u64>, u64>) -> EdgeListGraph> = 
    vec![
        directed_q::new::single_thread::get_q_digraph, 
        directed_q::new::mutex::get_q_digraph,
        directed_q::new::split_and_merge::get_q_digraph
        ];
    let names:Vec<String> = vec!["single_thread".to_string(), "mutex".to_string(), "split_and_merge".to_string()];
    let old_functions: Vec<for<'a, 'b> fn(usize, usize, usize, &'a Vec<Vec<Vec<u64>>>, &'b HashMap<Vec<u64>, u64>) -> EdgeListGraph> = 
    vec![
        directed_q::old::single_thread::get_q_digraph, 
        directed_q::old::mutex::get_q_digraph,
        directed_q::old::split_and_merge::get_q_digraph
        ];
    println!("==========================NEW======================");
    for (f, name) in zip(new_functions, names.clone())
    {
        let now = Instant::now();
        f( q, i, j, &flag_complex, &simplex_map);
        let elapsed = now.elapsed();
        println!("{} : {:.2?}", name, elapsed);
    }

    let now = Instant::now();
    directed_q::new::bottom_up::get_q_digraph(& g, q, i, j);
    let elapsed = now.elapsed();
    println!("Bottom Up : {:.2?}", elapsed);
    if mode == "full"
    {
        let now = Instant::now();
        directed_q::new::naive::get_q_digraph(g.clone(), q, i, j, &flag_complex);
        let elapsed = now.elapsed();
        println!("Naive: {:.2?}", elapsed);
    }
    println!("==========================OLD======================");
    for (f, name) in zip(old_functions, names)
    {
        let now = Instant::now();
        f( q, i, j, &flag_complex, &simplex_map);
        let elapsed = now.elapsed();
        println!("{} : {:.2?}", name, elapsed);
    }
    let now = Instant::now();
    directed_q::old::bottom_up::get_q_digraph(& g, q, i, j);    
    let elapsed = now.elapsed();
    println!("Bottom Up : {:.2?}", elapsed);
    if mode == "full"
    {
        let now = Instant::now();
        directed_q::old::naive::get_q_digraph(g.clone(), q, i, j, &flag_complex);
        let elapsed = now.elapsed();
        println!("Naive: {:.2?}", elapsed);
    }
}





