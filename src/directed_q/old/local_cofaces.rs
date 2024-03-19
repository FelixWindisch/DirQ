use array_tool::vec::Uniq;
use rayon::prelude::*;
use std::sync::Mutex;
use std::collections::HashMap;
use crate::graph::{EdgeMapGraph, DirectedGraphNew, EdgeListGraph};
use crate::graph::*;
use crate::util;
use itertools::Itertools;


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
        let simplex_inclusion = util::get_super_simplices(g, &q_simplex, simplex_map, &simplicial_family);
 
        for s in &simplex_inclusion
        {
            for t in &simplex_inclusion
            {
                //if *s == *t
                //{
                //    continue;
                //}
                let sigma = *s;
                let tau = *t;
                
                let sigma_i_cofaces:Vec<Node> = util::local_cofaces(g, simplicial_family.get(sigma as usize).unwrap(),  i as usize, simplex_map);
                let tau_j_cofaces:Vec<Node> = util::local_cofaces(g, simplicial_family.get(tau as usize).unwrap(),  j as usize, simplex_map);
        
                let mut q_graph_mut = q_graph.lock().unwrap();
                //q_graph_mut.add_edge(sigma, tau);    
                for from in sigma_i_cofaces
                {
                    for to in &tau_j_cofaces
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
    EdgeListGraph::new_disconnected(0)
    //q_graph
}