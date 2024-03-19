use rayon::prelude::*;
use std::sync::Mutex;
use std::collections::HashMap;
use crate::graph::{EdgeMapGraph, DirectedGraphNew, EdgeListGraph};
use crate::graph::*;
use crate::util;



pub fn get_q_digraph
(
    g: &EdgeMapGraph, 
    q:usize, i:usize, j:usize, 
    flag_complex:&Vec<Vec<Simplex>>,
    simplex_map: &HashMap<Simplex, Node>,
    
) -> EdgeListGraph 
{
    let simplicial_family: Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(|simplex: &Vec<Node>|{simplex.len() > q as usize}).collect();

    let q_graph = Mutex::new(EdgeListGraph::new_disconnected(0));

    let _enumerate = |q_simplex: Vec<Node>| {
        let simplex_i_cofaces:Vec<Node> = util::local_cofaces(g, &q_simplex,  i as usize, simplex_map);
        let simplex_j_cofaces:Vec<Node> = util::local_cofaces(g, &q_simplex,  j as usize, simplex_map);
        for i_coface in simplex_i_cofaces
        {
            for j_coface in &simplex_j_cofaces
            {
                let sigma = i_coface;
                let tau = *j_coface;
                
                let sigma_inclusion = util::get_super_simplices(g, &simplicial_family[sigma as usize], simplex_map, &simplicial_family);
                let tau_inclusion = util::get_super_simplices(g, &simplicial_family[tau as usize], simplex_map, &simplicial_family);
 
                let mut q_graph_mut = q_graph.lock().unwrap();
                q_graph_mut.add_edge(sigma, tau);    
                for from in sigma_inclusion
                {
                    for to in &tau_inclusion
                    {
                        if from != *to
                        {
                            q_graph_mut.add_edge(from, *to);
                        }
                    }
                }
                drop(q_graph_mut);
            }
        }
    };
    
    flag_complex[0].clone().into_par_iter().for_each(_enumerate);
    let q_graph = q_graph.lock().unwrap().to_owned();
    dbg!(q_graph.edges().len());

    q_graph
}