use array_tool::vec::Intersect;
extern crate array_tool;
use std::collections::{HashMap, HashSet};
use crate::graph::{EdgeMapGraph, Node};
pub use crate::graph::{Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;

/// Returns true if sigma is a face of tau
pub fn inclusion(sigma: &Simplex, tau: &Simplex)->bool
{
    //TODO: improve performance (merge-type algorithm?)
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

/// Modifies an integer vector by applying binary AND operation with another integer vector
fn binary_and(a:&mut Vec<u64>, b:&[u64])
{
    for (index, value) in b.iter().enumerate()
    {
        a[index] &= *value;
    }
}

/// Implementation of local coboundary search from Flagser.
/// On input of simplex σ and index i, returns vertices v in the original graph such that
/// (σ_0, ..., σ_{i-1}, v, σ_{i+1}, ..., σ_n) is in flag_complex
pub fn local_coboundaries<G:DirectedGraph + DirectedGraphNew + AdjacencyMatrixGraph>
(
    g: &G, 
    simplex: &Vec<Node>, 
    index: usize
) -> Vec<Node> 
{
    let dimension = simplex.len();
    let mut co_b: Vec<Node> = Vec::new();
    let mut c: Vec<crate::graph::Chunk> = vec![u64::MAX;g.row_len()];
    for j in 0..index
    {
        let b = g.edges_from_as_bitmap(simplex[j]);
        binary_and(&mut c, b);
    }

    for j in index..dimension {
        let b = g.edges_to_as_bitmap(simplex[j]);
        binary_and(&mut c, b);
    }

    for v in 0..g.nnodes() {
        let i = v / 64;
        let bit_offset = 1 << (v % 64);
        if c[i] & (bit_offset) != 0 {
            co_b.push(v as Node)
        }
    }
    return co_b;
}



/// On input of simplex and index, returns cof_index(simplex) as a vector of simplex IDs.
pub fn local_cofaces<G:DirectedGraph + DirectedGraphNew + AdjacencyMatrixGraph>
(
    g: &G, simplex: &Simplex, index: usize, 
    simplex_map: &HashMap<Simplex, Node>
) -> Vec<Node> 
{
    let co_b: Vec<Node> = local_coboundaries(g, simplex, index);
    let mut result: Vec<Node> = Vec::new();
    for co_boundary_vertex in co_b
    {
        let mut clone = simplex.clone();
        clone.insert(index, co_boundary_vertex);
        result.push(*simplex_map.get(&clone).unwrap());
    }
    return result;
}



/// On input of simplex and index, returns cof_index(simplex) as a vector of simplex IDs.
pub fn local_cofaces_full<G:DirectedGraph + DirectedGraphNew + AdjacencyMatrixGraph>
(
    g: &G, simplex: &Simplex, index: usize
) -> Vec<Simplex> 
{
    let co_b: Vec<Node> = local_coboundaries(g, simplex, index);
    let mut result: Vec<Simplex> = Vec::new();
    for co_boundary_vertex in co_b
    {
        let mut clone = simplex.clone();
        clone.insert(index, co_boundary_vertex);
        result.push(clone);
    }
    return result;
}

/// Returns the set of simplices that contain simplex as a vector of IDs
pub fn get_super_simplices(g: &EdgeMapGraph, simplex: &Vec<Node>, simplex_map: &HashMap<Simplex, Node>, simplicial_family: &Vec<Simplex>)-> Vec<Node> 
{
    let mut result:HashSet<Node> = HashSet::new();
    result.insert(*simplex_map.get(simplex).unwrap());
    // Rust closures cannot do recursion, so this is slightly ugly
    fn recurse(simplex: &Simplex, result: &mut HashSet<Node>, g: &EdgeMapGraph, simplicial_family: &Vec<Simplex>, simplex_map: &HashMap<Simplex, Node>) 
    {
        for index in 0..simplex.len()+1
        {

            let coface_ids = local_cofaces(g, &simplex,  index, simplex_map);
            result.extend(coface_ids.clone());
            for coface_id in coface_ids
            {
                let coface = &simplicial_family[coface_id as usize];
                recurse(&coface, result, g, simplicial_family, simplex_map);
            }
        }    
    }
    recurse(simplex, &mut result, g, simplicial_family, simplex_map);
    result.into_iter().collect()
}


/// Returns the set of simplices that contain simplex as a vector of IDs
pub fn get_super_simplices_full(g: &EdgeMapGraph, simplex: &Vec<Node>)-> Vec<Simplex> 
{
    let mut result:HashSet<Simplex> = HashSet::new();
    result.insert(simplex.clone());
    // Rust closures cannot do recursion, so this is slightly ugly
    fn recurse(simplex: &Simplex, result: &mut HashSet<Simplex>, g: &EdgeMapGraph) 
    {
        for index in 0..simplex.len()+1
        {

            let cofaces = local_cofaces_full(g, &simplex,  index);
            result.extend(cofaces.clone());
            for coface in cofaces
            {
                //let coface = &simplicial_family[coface_id as usize];
                recurse(&coface, result, g);
            }
        }    
    }
    recurse(simplex, &mut result, g);
    result.into_iter().collect()
}

pub fn get_super_simplices_with_inclusion(g: &EdgeMapGraph, simplex: &Vec<Node>)-> (Vec<Simplex>, HashSet<(Simplex, Simplex)>) 
{
    let mut inclusions:HashSet<(Simplex, Simplex)> = HashSet::new();
    let mut result:HashSet<Simplex> = HashSet::new();
    result.insert(simplex.clone());
    // Rust closures cannot do recursion, so this is slightly ugly
    fn recurse(simplex: &Simplex, result: &mut HashSet<Simplex>, inclusions: &mut HashSet<(Simplex, Simplex)>,  mut path: Vec<Simplex>, g: &EdgeMapGraph) 
    {
        path.push((*simplex).clone());
        for index in 0..simplex.len()+1
        {

            let cofaces = local_cofaces_full(g, &simplex,  index);
            result.extend(cofaces.clone());
            for coface in cofaces
            {

                for previous in path.clone()
                {
                    inclusions.insert((previous, coface.clone()));
                } 
                
                //let coface = &simplicial_family[coface_id as usize];
                recurse(&coface, result, inclusions, path.clone(), g);
            }
        }    
    }
    recurse(simplex, &mut result, &mut inclusions, Vec::new(), g);
    (result.into_iter().collect(), inclusions)
}






