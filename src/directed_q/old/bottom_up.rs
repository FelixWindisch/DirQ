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
            
            let (super_simplices, inclusions) = util::get_super_simplices_with_inclusion(g, &q_simplex);
            
            let mut q_graph_mut = result.lock().unwrap();
            q_graph_mut.extend(inclusions);
            drop(q_graph_mut);


            for sigma in super_simplices.clone()
            {
                for tau in super_simplices.clone()
                {
            
                    let sigma_i_cofaces:Vec<Simplex> = util::local_cofaces_full(g, &sigma,  i as usize);
                    let tau_j_cofaces:Vec<Simplex> = util::local_cofaces_full(g, &tau,  j as usize);
        
                    
                    let mut q_graph_mut = result.lock().unwrap();
                    for from in &sigma_i_cofaces
                    {
                        for to in &tau_j_cofaces
                        {
                            if *from != *to
                            {
                                q_graph_mut.insert(((*from).clone(), (*to).clone()));
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