use crate::graph::*;

use rayon::prelude::*;
use std::sync::Mutex;
use std::collections::HashMap;
use crate::graph::{Simplex, Edge, Node, DirectedGraphNew, EdgeListGraph};


use itertools::Itertools;

/// Implementation of efficient algorithm using Mutual Exclusion Parallelization to compute the 
/// (q, i, j)-digraph according to the new definition
pub fn get_q_digraph( q:usize, i:usize, j:usize, flag_complex: &Vec<Vec<Simplex>>, simplex_map: &HashMap<Simplex, Node>) -> EdgeListGraph
{
    // Cofaces and inclusions can be computed in parallel
    let (i_j_cofaces, (q_plus_one_inclusions, inclusion_edges))
    = rayon::join(
        ||{get_i_j_cofaces( q, i, j,  flag_complex[0].len(), &flag_complex[1][..], &simplex_map)},
        ||{get_inclusion_edges( q, i, j, flag_complex[1].len(), &flag_complex, &simplex_map)});

    let mut q_graph = get_q_near_graph( q, i, j,   flag_complex[0].len(), &i_j_cofaces, &q_plus_one_inclusions);
    for edge in inclusion_edges
    {
        q_graph.add_edge(edge[0], edge[1]);
    }
    dbg!(q_graph.edges().len());
    q_graph
}



pub fn get_i_j_cofaces
(
    _q:usize, i:usize, j:usize,  
    number_q_simplices:usize, 
    q_plus_one_simplices: &[Vec<Node>], 
    simplex_map: &HashMap<Simplex, Node>
) -> Vec<(Vec<Node>, Vec<Node>)> 
{
    let result = Mutex::new(vec![(Vec::new(), Vec::new()); number_q_simplices]);
    let enumerate = |simplex: &Vec<Node>| {

            let mut i_boundary =simplex.clone();
            i_boundary.remove(i as usize);
            let mut j_boundary = simplex.clone();
            j_boundary.remove(j as usize);
            
            let i_q_simplex = *simplex_map.get(&i_boundary).unwrap();
            let j_q_simplex = *simplex_map.get(&j_boundary).unwrap();
            let index = *simplex_map.get(simplex).unwrap();
            // Code block here is important, so the mutex guard goes out of scope
            {
                let mut result_guard = result.lock().unwrap();
                result_guard[i_q_simplex as usize].0.push(index);
                result_guard[j_q_simplex as usize].1.push(index);
            }

    };
    q_plus_one_simplices.into_par_iter().for_each(enumerate);
    let result = result.lock().unwrap().to_owned();
    result
}



pub fn get_inclusion_edges
(
    q:usize, _i:usize, _j:usize, 
    number_q_plus_one_simplices:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> (Vec<Vec<Node>>, Vec<Edge>) 
{
    let q_plus_one_inclusions = Mutex::new(vec![Vec::new(); number_q_plus_one_simplices]);
    let inclusion_edges = Mutex::new(Vec::new());
    let enumerate = |simplex: &Vec<Node>| {
        let simplex_id = *simplex_map.get(simplex).unwrap();
        for face_len in q as usize..simplex.len()
        {
            let faces : Vec<Vec<Node>> = simplex.clone().into_iter().combinations(face_len).collect();
            for face in faces
            {
                //TODO: This is silly. Roll your own combinations() for Node?
                let face: Vec<Node> = face.into_iter().map(|x| x).collect();
                
                if face.len() > q as usize 
                {
                    let face_id = *simplex_map.get(&face).unwrap();
                    if face.len() == (q+2) as usize
                    {
                        let mut q_plus_one_simplex_inclusions = q_plus_one_inclusions.lock().unwrap();
                        q_plus_one_simplex_inclusions[(face_id - flag_complex[0].len() as Node) as usize].push(simplex_id);
                    }
                    let mut inc_edges_mut = inclusion_edges.lock().unwrap();
                    inc_edges_mut.push([face_id as Node, simplex_id as Node]);
                }
            }
        }
    };
    flag_complex.into_par_iter().skip(1).flatten().for_each(enumerate);
    let mut q_plus_one_inclusions = q_plus_one_inclusions.lock().unwrap().to_owned();
    let inclusion_edges = inclusion_edges.lock().unwrap().to_owned();
    for i in 0..number_q_plus_one_simplices
    {
        q_plus_one_inclusions[i].push(i as Node + flag_complex[0].len() as Node);
    }
    (q_plus_one_inclusions, inclusion_edges)
}



pub fn get_q_near_graph
(
    _q:usize, _i:usize, _j:usize, 
    number_q_simplices:usize, 
    i_j_cofaces: &Vec<(Vec<Node>, Vec<Node>)>, 
    q_plus_one_simplex_inclusions : &Vec<Vec<Node>>
) -> EdgeListGraph 
{
    let q_graph = Mutex::new(EdgeListGraph::new_disconnected(0));

    let enumerate = |q_simplex_id: usize| {
        let simplex_i_cofaces:&Vec<Node> = &i_j_cofaces[q_simplex_id as usize].0;
        let simplex_j_cofaces:&Vec<Node> = &i_j_cofaces[q_simplex_id as usize].1;

        for i_coboundary in simplex_i_cofaces
        {
            for j_coboundary in simplex_j_cofaces
            {

                let sigma = *i_coboundary;
                let tau = *j_coboundary;
                
                let tau_inclusion = &q_plus_one_simplex_inclusions[(tau - number_q_simplices as Node) as usize];
                let sigma_inclusion = &q_plus_one_simplex_inclusions[(sigma  - number_q_simplices as Node) as usize];
                
                let mut q_graph_guard = q_graph.lock().unwrap();
                q_graph_guard.add_edge(sigma, tau);    
                for from in sigma_inclusion
                {
                    for to in tau_inclusion
                    {
                        if *from != *to
                        {
                            q_graph_guard.add_edge(*from, *to);
                        }
                    }
                }
                drop(q_graph_guard);
            }
        }
    };
    
    (0..number_q_simplices).into_par_iter().for_each(enumerate);
    let q_graph = q_graph.lock().unwrap().to_owned();
    q_graph
}
