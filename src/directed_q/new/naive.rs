extern crate array_tool;


use crate::graph::EdgeListGraph;
use array_tool::vec::Intersect;


use crate::graph::EdgeMapGraph;

pub use crate::graph::{Node, Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;
use itertools::Itertools;
use crate::util;





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
            else if simplicial_family[sigma].len() > q as usize && simplicial_family[tau].len() > q as usize && are_q_near_new(simplicial_family[sigma].clone(),simplicial_family[tau].clone(), q, i, j)
            {
                q_graph.add_edge(sigma as Node,tau as Node);
            }
        }
    }
    dbg!(q_graph.edges().len());
    q_graph
}


fn are_q_near_new(sigma: Vec<Node>, tau:Vec<Node>, q:usize, i:usize, j:usize) ->  bool
{
    let common_vertices = sigma.intersect(tau.clone());
    if common_vertices.len() > q as usize
    {
        let sigma_faces = get_q_faces_coboundary(sigma.clone(), q, i);
        let tau_faces = get_q_faces_coboundary(tau.clone(), q, j);
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

// calculates cob_i(simplex)
fn get_q_faces_coboundary(simplex: Vec<Node>, q:usize, i_hat:usize) -> Vec<Vec<Node>>
{
    let faces:Vec<Vec<Node>> = simplex.into_iter().combinations((q+2) as usize).collect();
    let mut result:Vec<Vec<Node>> = Vec::new();
    for mut face in faces
    {
        face.remove(i_hat as usize);
        result.push(face);
    }
    return  result;
}