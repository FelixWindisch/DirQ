

extern crate array_tool;
use std::collections::HashMap;

use rayon::prelude::*;
use std::sync::Mutex;

use crate::graph::EdgeListGraph;



use crate::graph::{Edge, Node, Simplex};

pub use crate::graph::{DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;
use itertools::Itertools;



/// Implementation of efficient algorithm using Mutual Exclusion parallelization to compute the 
/// (q, i, j)-digraph according to the original definition
pub fn get_q_digraph
(
    q:usize, i:usize, j:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> EdgeListGraph
{

    // We can find the cofaces and inclusions in parallel
    let (i_j_cofaces, (q_simplex_inclusions, inclusion_edges))
    = rayon::join(
        ||{get_all_i_j_cofaces( q, i, j, &flag_complex, simplex_map)},
        ||{get_inclusion_edges( q, i, j, &flag_complex, &simplex_map)});
    get_q_near_graph( q, i, j,  i_j_cofaces, inclusion_edges, q_simplex_inclusions)
}


pub fn get_all_i_j_cofaces
(
    q:usize, i:usize, j:usize,  
    _flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> Vec<(Vec<Node>, Vec<Node>)> 
{
    let number_simplices = simplex_map.len();
    let result = Mutex::new(vec![(Vec::new(), Vec::new()); number_simplices]);
    let enumerate = |simplex: &Vec<Node>| {
        if simplex.len() > q + 1
        {
            let mut i_boundary :Vec<Node> = simplex.clone();
            i_boundary.remove(i as usize);
            let mut j_boundary :Vec<Node> = simplex.clone();
            j_boundary.remove(j as usize);
            
            let i_face = *simplex_map.get(&i_boundary).unwrap();
            let j_face = *simplex_map.get(&j_boundary).unwrap();
            let index = *simplex_map.get(simplex).unwrap();
            {
                let mut result_guard = result.lock().unwrap();
                result_guard[i_face as usize].0.push(index);
                result_guard[j_face as usize].1.push(index);
            }
        }
    };
    simplex_map.keys().collect::<Vec<_>>().into_par_iter().for_each(enumerate);
    let result = result.lock().unwrap().to_owned();
    result
}

pub fn get_inclusion_edges
(
    q:usize, _i:usize, _j:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> (Vec<Vec<Node>>, Vec<Edge>) 
{
    let q_inclusions = Mutex::new(vec![Vec::new(); flag_complex[0].len()]);
    let inclusion_edges = Mutex::new(Vec::new());
    let enumerate = |simplex: &Vec<Node>| {
        let simplex_id = *simplex_map.get(simplex).unwrap();
        for face_len in q as usize..simplex.len()
        {
            let faces : Vec<Vec<Node>> = simplex.clone().into_iter().combinations(face_len).collect();
            for face in faces
            {
                //This is silly. Roll your own combinations() for Node?
                let face: Vec<Node> = face.into_iter().map(|x| x).collect();
                
                if face.len() > q as usize 
                {
                    let face_id = *simplex_map.get(&face).unwrap();
                    if face.len() == (q+1) as usize
                    {
                        let mut q_simplex_inclusions = q_inclusions.lock().unwrap();
                        q_simplex_inclusions[face_id as usize].push(simplex_id);
                    }
                    let mut inc_edges_mut = inclusion_edges.lock().unwrap();
                    inc_edges_mut.push([face_id as Node, simplex_id as Node]);
                }
            }
        }
    };
    flag_complex.into_par_iter().skip(1).flatten().for_each(enumerate);
    let mut q_inclusions = q_inclusions.lock().unwrap().to_owned();
    let inclusion_edges = inclusion_edges.lock().unwrap().to_owned();
    for i in 0..flag_complex[0].len()
    {
        q_inclusions[i].push(i as Node);
    }
    (q_inclusions, inclusion_edges)
}


pub fn get_q_near_graph
(
    _q:usize, _i:usize, _j:usize, 
    i_j_cofaces: Vec<(Vec<Node>, Vec<Node>)>, 
    inclusion_edges: Vec<Edge>, 
    q_simplex_inclusions : Vec<Vec<Node>>
) -> EdgeListGraph
{   
    let q_graph = Mutex::new(EdgeListGraph::new_disconnected(0));
    let enumerate = |q_simplex_id: usize| 
    {
        let q_simplex_inclusion:&Vec<Node> = &q_simplex_inclusions[q_simplex_id as usize];
        for sigma in q_simplex_inclusion
        {
            for tau in q_simplex_inclusion
            {
                // self-loops are disregarded
                //if sigma == tau
                //{
                //    continue;
                //}
                let sigma_cofaces = &i_j_cofaces[*sigma as usize].0;
                let tau_cofaces = &i_j_cofaces[*tau as usize].1;

                let mut q_graph_guard = q_graph.lock().unwrap();

                for from in sigma_cofaces
                {
                    for to in tau_cofaces
                    {
                        if *from != *to
                        {
                            q_graph_guard.add_edge(*from, *to);
                        }
                    }
                }
            }
        }
    };
    (0..q_simplex_inclusions.len()).into_par_iter().for_each(enumerate);
    let mut q_near_graph = q_graph.lock().unwrap().to_owned();
    for edge in &inclusion_edges
    {
        q_near_graph.add_edge(edge[0], edge[1]);
    }
    dbg!(q_near_graph.edges().len());
    q_near_graph
}




