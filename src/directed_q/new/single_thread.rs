extern crate array_tool;
use crate::DirectedGraph;
use std::collections::HashMap;
use crate::graph::{Node, Simplex, Edge, DirectedGraphNew, EdgeListGraph};
use itertools::Itertools;

/// Single threaded implementation of efficient algorithm to compute the 
/// (q, i, j)-digraph according to the new definition
pub fn get_q_digraph
(
    q:usize, i:usize, j:usize, 
    flag_complex: &Vec<Vec<Simplex>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> EdgeListGraph
{
    let (q_plus_one_super_simplices, inclusion_edges) = get_inclusion_edges( q, i, j,   &flag_complex, &simplex_map);
    let i_j_cofaces = get_i_j_cofaces( q, i, j,  flag_complex[0].len(), flag_complex[1].clone(), &simplex_map);

    let q_graph:EdgeListGraph = get_q_near_graph(  q, i, j,  i_j_cofaces, inclusion_edges, q_plus_one_super_simplices);

    dbg!(q_graph.edges().len());
    q_graph
}


pub fn get_i_j_cofaces
(
    _q:usize, i:usize, j:usize, 
    number_q_simplices:usize, 
    q_plus_one_simplices: Vec<Vec<Node>>, 
    simplex_map: &HashMap<Simplex, Node>
) -> Vec<(Vec<Node>, Vec<Node>)> 
{
    let mut result = vec![(Vec::new(), Vec::new()); number_q_simplices];
    for simplex in q_plus_one_simplices{

        let mut i_boundary :Vec<Node> = simplex.clone();
        i_boundary.remove(i as usize);
        let mut j_boundary :Vec<Node> = simplex.clone();
        j_boundary.remove(j as usize);
        
        let i_q_simplex = *simplex_map.get(&i_boundary).unwrap();
        let j_q_simplex = *simplex_map.get(&j_boundary).unwrap();
        let index = *simplex_map.get(&simplex).unwrap();
        
        result[i_q_simplex as usize].0.push(index);
        result[j_q_simplex as usize].1.push(index);
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
    let mut inclusion_graph = Vec::<Edge>::new();
    // inclusion edges
    let mut q_plus_one_simplex_inclusions:Vec<Vec<Node>> = vec![Vec::new(); flag_complex[1].len()];
    for q_plus_one_simplex_id in 0..flag_complex[1].len()
    {
        q_plus_one_simplex_inclusions[q_plus_one_simplex_id as usize].push(q_plus_one_simplex_id as Node + flag_complex[0].len() as Node)
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
                        if face.len() == (q+2) as usize
                        {
                            q_plus_one_simplex_inclusions[(face_id - flag_complex[0].len() as Node) as usize].push(simplex_id);
                        }
                        inclusion_graph.push([face_id as Node, simplex_id as Node]);
                    }
                }
            }
        }
    }
    return (q_plus_one_simplex_inclusions, inclusion_graph)
}

pub fn get_q_near_graph
(
    _q:usize, _i:usize, _j:usize, 
    i_j_cofaces: Vec<(Vec<Node>, Vec<Node>)>, 
    inclusion_edges: Vec<Edge>, 
    q_plus_one_simplex_inclusions : Vec<Vec<Node>>
) -> EdgeListGraph
{   
    let mut q_near_graph = EdgeListGraph::new_disconnected(0);
    let number_q_simplices = i_j_cofaces.len();
    for q_simplex_id in 0..number_q_simplices
    {
        let simplex_i_cofaces:&Vec<Node> = &i_j_cofaces[q_simplex_id as usize].0;
        let simplex_j_cofaces:&Vec<Node> = &i_j_cofaces[q_simplex_id as usize].1;

        for i_coboundary in simplex_i_cofaces
        {
            for j_coboundary in simplex_j_cofaces
            {
                let sigma = *i_coboundary;
                let tau = *j_coboundary;
                q_near_graph.add_edge(sigma, tau);
                
                let tau_inclusion = &q_plus_one_simplex_inclusions[(tau - number_q_simplices as Node) as usize];
                let sigma_inclusion = &q_plus_one_simplex_inclusions[(sigma - number_q_simplices as Node) as usize];
                for from in sigma_inclusion
                {
                    for to in tau_inclusion
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
    q_near_graph
}