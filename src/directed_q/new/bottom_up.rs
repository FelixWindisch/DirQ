use rayon::prelude::*;
use std::sync::Mutex;
use std::collections::{HashMap, HashSet};
use crate::graph::{EdgeMapGraph, DirectedGraphNew, EdgeListGraph};
use crate::graph::*;
use crate::util;
use std::cmp::max;


pub fn get_q_digraph(g: &EdgeMapGraph, q:usize, i:usize, j:usize) -> HashSet<(Simplex, Simplex)> {
    let mut result = Mutex::new(HashSet::new());
    let enumerate = |_i: &mut u64, simplex: &[Node]| 
    {
        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(simplex);
        if x.len() == (q+1)
        {
            
            let q_simplex = x;
            
            let (_, inclusions) = util::get_super_simplices_with_inclusion(g, &q_simplex);
            
            let mut q_graph_mut = result.lock().unwrap();
            //q_graph_mut.extend(inclusions);
            drop(q_graph_mut);


            let simplex_i_cofaces:Vec<Simplex> = util::local_cofaces_full(g, &q_simplex,  i as usize);
            let simplex_j_cofaces:Vec<Simplex> = util::local_cofaces_full(g, &q_simplex,  j as usize);
            for i_coface in simplex_i_cofaces
            {
            for j_coface in simplex_j_cofaces.clone()
                {
                    let sigma = i_coface.clone();
                    let tau = j_coface;
            
                    let sigma_inclusion = util::get_super_simplices_full(g, &sigma);
                    let tau_inclusion = util::get_super_simplices_full(g, &tau);
                    
                    let mut q_graph_mut = result.lock().unwrap();
                    //q_graph_mut.insert((sigma, tau));    
                    for from in &sigma_inclusion
                    {
                        for to in &tau_inclusion
                        {
                            if *from != *to
                            {
                                //q_graph_mut.insert(((*from).clone(), (*to).clone()));
                            }
                        }
                    }
                    drop(q_graph_mut);
                }
            }

        }
    };
    
    crate::complex::for_each_cell_par(g, &enumerate, 0, g.nnodes()).reduce(||{0}, &|x, y|{x + y});
    let q_graph = result.lock().unwrap().to_owned();
    dbg!(q_graph.len());
    q_graph
}