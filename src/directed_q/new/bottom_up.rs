use std::sync::Mutex;
use std::collections::HashSet;
use rayon::iter::ParallelIterator;

use crate::graph::EdgeMapGraph;
use crate::graph::*;
use crate::util;


pub fn get_q_digraph(g: &EdgeMapGraph, q:usize, i:usize, j:usize) -> HashSet<(Simplex, Simplex)> {
    let result = Mutex::new(HashSet::new());
    let enumerate = |_i: &mut u64, simplex: &[Node]| 
    {
        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(simplex);
        {
            
            let q_simplex = x;
            
            let (_, inclusions) = util::get_super_simplices_with_inclusion(g, &q_simplex);
            
            let mut q_graph_mut = result.lock().unwrap();
            q_graph_mut.extend(inclusions);
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
                    q_graph_mut.insert((sigma, tau));    
                    for from in &sigma_inclusion
                    {
                        for to in &tau_inclusion
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
    
    let _y = crate::complex::for_each_cell_par(g, &enumerate, q, q).reduce(|| 0, |x,_y| x);
    let q_graph = result.lock().unwrap().to_owned();
    dbg!(q_graph.len());
    q_graph
}