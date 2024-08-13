pub mod new;
pub mod original;
use crate::complex;
use crate::graph;
use crate::graph_io;
use crate::directed_q;
pub use graph::{EdgeMapGraph, EdgeListGraph, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};

pub type Graph = graph::EdgeMapGraph;
use std::time::Instant;
use std::collections::HashMap;
use core::iter::zip;
use clap::Parser;
use std::path::Path;


#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
#[derive(Default)]
pub struct Args{
    // Configuration
    /// flag input file location
    #[clap(long, default_value="test.flag")]
    pub input_file: String,

    #[clap(long, default_value="out/")]
    pub output_folder: String,


    /// direction_from
    #[clap(short, long, default_value="0")]
    pub i: usize,

    /// direction_to
    #[clap(short, long, default_value="2")]
    pub j: usize,

    /// q
    #[clap(short, long, default_value="2")]
    pub q: usize,

    /// maximal simplex dimension
    #[clap(short, long, default_value="10000000")]
    pub max_dimension: usize,

    /// definition: original (Henries) or new (Florians/Felix')
    #[clap(short, long, possible_values=&["original", "new"], default_value="new")]
    pub definition: String,


    /// parallelization method
    #[clap(short, long, possible_values=&["none", "mutex", "split_and_merge", "local_cofaces", "bottom_up", "naive"], default_value="none")]
    pub parallelization: String,

    #[clap(short, long, possible_values=&["false", "parallel", "true", "full"], default_value="false")]
    pub benchmark: String
}
pub fn run(args:Args)
{
    let i = args.i;
    let j = args.j;
    let q = args.q;
    let max_dimension = args.max_dimension;
    assert!(max_dimension > q);
    let par = args.parallelization;
    let definition = args.definition;
    let input_file  = args.input_file;
    let output_folder = args.output_folder;
    let bm = args.benchmark;

    let path = Path::new(&input_file);
    let filename = path.file_name().unwrap().to_str().unwrap();
    let mut export_file = output_folder; 
    export_file.push_str("/");
    export_file.push_str(filename.split('.').next().unwrap());
    export_file.push_str( &(format!("_q{}i{}j{}_{}", args.q, args.i, args.j, definition)));
    if par == "bottom_up"
    {
        export_file.push_str(".q_graph");
    }
    else
    {
        export_file.push_str(".flag");
    }
    let g: EdgeMapGraph = if input_file == "random.flag"
    {
        EdgeMapGraph::gen_seo_er(278, 0.028491286 *2.0, &mut rand::thread_rng())
    }
    else
    {
         graph_io::read_flag_file(&input_file)
    };
    if bm == "true" || bm == "full" || bm == "parallel"
    {
        benchmark(g, q, i, j, bm);
        return;
    } 
    if par == "bottom_up"
    {
        if definition == "original"
        {
            let result = directed_q::original::bottom_up::get_q_digraph(&g, q, i, j);
            graph_io::save_flag_file_edge_list(&export_file, result);
        }
        else
        {
            let result = directed_q::new::bottom_up::get_q_digraph(&g, q, i, j);
            graph_io::save_flag_file_edge_list(&export_file, result);
            
        }
            
        
        return;
    }

    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q, max_dimension);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);
    assert!(flag_complex.len() > 1);
    if flag_complex[1].len() == 0
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
                //"bottom_up" => {directed_q::new::bottom_up::get_q_digraph(&g, q, i, j)}
                _ => {directed_q::new::single_thread::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
            }
        }
        "original" =>
        {
            match par.as_str()
            {
                "naive" => {directed_q::original::naive::get_q_digraph(g, q, i, j, &flag_complex)}
                "mutex" => {directed_q::original::mutex::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "split_and_merge" => {directed_q::original::split_and_merge::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
                "local_cofaces" => {directed_q::original::local_cofaces::get_q_digraph( &g, q, i, j, &flag_complex, &simplex_map)}
                //"bottom_up" => {directed_q::original::bottom_up::get_q_digraph(&g, q, i, j)}
                _ => {directed_q::original::single_thread::get_q_digraph( q, i, j, &flag_complex,  &simplex_map)}
            }
        }
        _ => {directed_q::new::single_thread::get_q_digraph( q, i, j, &flag_complex, &simplex_map)}
    };
    
    
    graph_io::save_flag_file(&export_file, &q);
}


fn benchmark(g: EdgeMapGraph, q:usize, i:usize, j:usize, mode: String)
{
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q, 100000000000);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);

    let new_functions: Vec<for<'a, 'b> fn(usize, usize, usize, &'a Vec<Vec<Vec<u64>>>, &'b HashMap<Vec<u64>, u64>) -> EdgeListGraph> = 
    vec![
        directed_q::new::single_thread::get_q_digraph, 
        directed_q::new::mutex::get_q_digraph,
        directed_q::new::split_and_merge::get_q_digraph
        ];
    let names:Vec<String> = vec![ "single_thread".to_string(), "mutex".to_string(), "split_and_merge".to_string()];
    let original_functions: Vec<for<'a, 'b> fn(usize, usize, usize, &'a Vec<Vec<Vec<u64>>>, &'b HashMap<Vec<u64>, u64>) -> EdgeListGraph> = 
    vec![
        directed_q::original::single_thread::get_q_digraph, 
        directed_q::original::mutex::get_q_digraph,
        directed_q::original::split_and_merge::get_q_digraph
        ];
    println!("==========================NEW======================");
    for (f, name) in zip(new_functions, names.clone())
    {
        let now = Instant::now();
        f( q, i, j, &flag_complex, &simplex_map);
        let elapsed = now.elapsed();
        println!("{} : {:.2?}", name, elapsed);
    }
    if mode != "parallel"
    {
        let now = Instant::now();
        directed_q::new::bottom_up::get_q_digraph(& g, q, i, j);
        let elapsed = now.elapsed();
        println!("Bottom Up : {:.2?}", elapsed);
    }
    if mode == "full"
    {
        let now = Instant::now();
        directed_q::new::naive::get_q_digraph(g.clone(), q, i, j, &flag_complex);
        let elapsed = now.elapsed();
        println!("Naive: {:.2?}", elapsed);
    }
    println!("=======================ORIGINAL====================");
    for (f, name) in zip(original_functions, names)
    {
        let now = Instant::now();
        f( q, i, j, &flag_complex, &simplex_map);
        let elapsed = now.elapsed();
        println!("{} : {:.2?}", name, elapsed);
    }
    if mode != "parallel"
    {
        let now = Instant::now();
        directed_q::original::bottom_up::get_q_digraph(& g, q, i, j);    
        let elapsed = now.elapsed();
        println!("Bottom Up : {:.2?}", elapsed);
    }
    if mode == "full"
    {
        let now = Instant::now();
        directed_q::original::naive::get_q_digraph(g.clone(), q, i, j, &flag_complex);
        let elapsed = now.elapsed();
        println!("Naive: {:.2?}", elapsed);
    }
}
