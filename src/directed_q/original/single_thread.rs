extern crate array_tool;
use std::collections::HashMap;


use crate::graph::EdgeListGraph;



use crate::graph::{Edge, Node, Simplex};

pub use crate::graph::{DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;
use itertools::Itertools;



/// Single-threaded implementation of efficient algorithm to compute the 
/// (q, i, j)-digraph according to the original definition
pub fn get_q_digraph
(
    q:usize, i:usize, j:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> EdgeListGraph
{
    let i_j_cofaces = get_all_i_j_cofaces( q, i, j, &flag_complex, simplex_map);
    let (q_simplex_inclusions, inclusion_edges) = get_inclusion_edges( q, i, j, &flag_complex, &simplex_map);
    get_q_near_graph( q, i, j, simplex_map.len(), i_j_cofaces, inclusion_edges, q_simplex_inclusions)
}


pub fn get_all_i_j_cofaces
(
    q:usize, i:usize, j:usize,  
    _flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> Vec<(Vec<Node>, Vec<Node>)> 
{
    let number_simplices = simplex_map.len();
    let mut result = vec![(Vec::new(), Vec::new()); number_simplices];
    for simplex in simplex_map.keys()
    {
        if simplex.len() > q + 1
        {
            let mut i_boundary :Vec<Node> = simplex.clone();
            i_boundary.remove(i as usize);
            let mut j_boundary :Vec<Node> = simplex.clone();
            j_boundary.remove(j as usize);
            
            let i_face = *simplex_map.get(&i_boundary).unwrap();
            let j_face = *simplex_map.get(&j_boundary).unwrap();
            let index = *simplex_map.get(simplex).unwrap();
            
            result[i_face as usize].0.push(index);
            result[j_face as usize].1.push(index);
        }
    }
    result
}

pub fn get_inclusion_edges
(
    q:usize, _i:usize, _j:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> (Vec<Vec<Node>>, Vec<Edge>)
{
    let mut inclusion_edges = Vec::<Edge>::new();
    // inclusion edges
    let mut q_simplex_inclusions:Vec<Vec<Node>> = vec![Vec::new(); flag_complex[0].len()];
    for q_simplex_id in 0..flag_complex[0].len()
    {
        q_simplex_inclusions[q_simplex_id as usize].push(q_simplex_id as Node);
    }

    for dimension in flag_complex.iter().skip(1).rev()
    {
        for simplex in dimension
        {
            let simplex_id = *simplex_map.get(simplex).unwrap();
            for face_len in q as usize..simplex.len()
            {
                let faces : Vec<Vec<&Node>> = simplex.into_iter().combinations(face_len).collect();
                for face in faces
                {
                    //This is silly. Roll your own combinations() for Node?
                    let face: Vec<Node> = face.into_iter().map(|x| *x).collect();
                    
                    if face.len() > q as usize 
                    {
                        let face_id = *simplex_map.get(&face).unwrap();
                        if face.len() == (q+1) as usize
                        {
                            q_simplex_inclusions[(face_id as Node) as usize].push(simplex_id);
                        }
                        inclusion_edges.push([face_id as Node, simplex_id as Node]);
                    }
                }
            }
        }
    }
    return (q_simplex_inclusions, inclusion_edges)
}


pub fn get_q_near_graph
(
    _q:usize, _i:usize, _j:usize, 
    number_simplices: usize,
    i_j_cofaces: Vec<(Vec<Node>, Vec<Node>)>, 
    inclusion_edges: Vec<Edge>, 
    q_simplex_inclusions : Vec<Vec<Node>>
) -> EdgeListGraph
{   
    let mut q_near_graph = EdgeListGraph::new_disconnected(number_simplices as usize);
    for q_simplex_id in 0..q_simplex_inclusions.len()
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

                for from in sigma_cofaces
                {
                    for to in tau_cofaces
                    {
                        if *from != *to
                        {
                            q_near_graph.add_edge(*from, *to);
                        }
                    }
                }
            }
        }
    }
    for edge in &inclusion_edges
    {
        q_near_graph.add_edge(edge[0], edge[1]);
    }
    dbg!(q_near_graph.edges().len());
    q_near_graph
}