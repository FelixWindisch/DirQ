extern crate array_tool;
use crate::graph::EdgeListGraph;
use array_tool::vec::Intersect;
use crate::graph::EdgeMapGraph;
pub use crate::graph::{Node, Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;
use itertools::Itertools;
use crate::util;




/// Single-threaded implementation of naive algorithm to compute the 
/// (q, i, j)-digraph according to the original definition
pub fn inclusion(sigma: &Vec<Node>, tau: &Vec<Node>)->bool
{
    let common_vertices = sigma.intersect(tau.to_vec());
    if common_vertices.len() == sigma.len()
    {
        let mut sorted_tau = common_vertices.clone();
        sorted_tau.sort_by(|a, b| tau.iter().position(|&x| x == *a).partial_cmp(&tau.iter().position(|&x| x == *b)).unwrap() ); 
        let mut sorted_sigma = common_vertices.clone();
        sorted_sigma.sort_by(|a, b| sigma.iter().position(|&x| x == *a).partial_cmp(&sigma.iter().position(|&x| x == *b)).unwrap() ); 
        
        return sorted_tau == sorted_sigma;
    }
    false
}


pub fn get_q_digraph(_graph: EdgeMapGraph, q:usize, i:usize, j:usize, flag_complex: &Vec<Vec<Simplex>>) -> EdgeListGraph
{
    let f = |simplex: &Vec<Node>|{simplex.len() > q as usize};
    let simplicial_family:Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(f).collect();
    
    let mut q_graph = EdgeListGraph::new_disconnected(simplicial_family.len());
    for sigma in 0..simplicial_family.len()
    {
        for tau in 0..simplicial_family.len()
        {
            if sigma == tau
            {
                continue;
            }
            else if simplicial_family[sigma].len() < simplicial_family[tau].len() && util::inclusion(&simplicial_family[sigma],&simplicial_family[tau])
            {
                q_graph.add_edge(sigma as Node,tau as Node);
            }
            else if simplicial_family[sigma].len() > q as usize && simplicial_family[tau].len() > q as usize && are_q_near_old(q,simplicial_family[sigma].clone(),simplicial_family[tau].clone(),i, j)
            {
                q_graph.add_edge(sigma as Node,tau as Node);
            }
        }
    }
    dbg!(q_graph.edges().len());
    q_graph
}



fn are_q_near_old(q:usize, sigma: Vec<Node>, tau:Vec<Node>, i:usize, j:usize) ->  bool
{
    let common_vertices = sigma.intersect(tau.clone());


    if common_vertices.len() > q as usize
    {
        let mut sigma_ = sigma.clone();
        sigma_.remove(i as usize);
        let mut tau_ = tau.clone();
        tau_.remove(j as usize);
        let sigma_faces = get_q_faces(sigma_, q);
        let tau_faces = get_q_faces(tau_, q);
        for sigma_face in sigma_faces
        {
            if tau_faces.contains(&sigma_face)
            {
                return true;
            }
        }
    }
    return false;
}

// calculates δ_i(simplex)
fn get_q_faces(simplex: Vec<Node>, q:usize) -> Vec<Vec<Node>>
{
    let faces:Vec<Vec<Node>> = simplex.into_iter().combinations((q+1) as usize).collect();
    let mut result:Vec<Vec<Node>> = Vec::new();
    for face in faces
    {
        result.push(face);
    }
    return  result;
}