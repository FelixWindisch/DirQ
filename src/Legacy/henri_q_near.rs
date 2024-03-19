extern crate array_tool;
use std::collections::HashMap;

use crate::graph::EdgeListGraph;



use crate::graph::{Node, Edge, EdgeMapGraph};

pub use crate::graph::{Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;
use itertools::Itertools;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};





use rayon::iter::IndexedParallelIterator;



pub fn get_inclusion_with_q<G: DirectedGraph + Sync + ?Sized>(_graph: &G, _i:usize, _j:usize, q:usize, number_q_simplices:usize, flag_complex: &Vec<Vec<Simplex>>, simplex_map: &HashMap<Simplex, Node>) -> (Vec<Vec<Node>>, Vec<Edge>) {
    let enumerate = |mut result: (Vec<Vec<Node>>, Vec<Edge>), simplices: &Vec<Vec<Node>>| {
        for simplex in simplices
        {
            let simplex_id = *simplex_map.get(simplex).unwrap();
            for face_len in q as usize..simplex.len()
            {
                let faces : Vec<Vec<Node>> = simplex.clone().into_iter().combinations(face_len).collect();
                for face in faces
                {
                    if face.len() > q as usize 
                    {
                        let face_id = *simplex_map.get(&face).unwrap();
                        if face.len() == (q+1) as usize
                        {
                            result.0[face_id as usize].push(simplex_id);
                        }
                        result.1.push([face_id as Node, simplex_id as Node]);
                    }
                }
            }
        }
        result
    };
    
    let combine = |right:  (Vec<Vec<Node>>, Vec<Edge>), left :  (Vec<Vec<Node>>, Vec<Edge>)|
    {
        let mut bigger_edges;
        let smaller_edges;
        let mut bigger_inclusion;
        let smaller_inclusion;
        if left.1.len() > right.1.len()
        {
            (bigger_inclusion, bigger_edges) = left;
            (smaller_inclusion, smaller_edges) = right;
        }
        else
        {
            (bigger_inclusion, bigger_edges) = right;
            (smaller_inclusion, smaller_edges) = left;
        }
        let _index = 0;
        bigger_edges.extend(smaller_edges);
        for (index, inclusion) in smaller_inclusion.iter().enumerate()
        {
            bigger_inclusion[index].extend(inclusion);
        }
        (bigger_inclusion, bigger_edges)
    };

    let mut result = flag_complex.into_par_iter().skip(1).chunks(100).flatten().fold(|| (vec![Vec::new(); number_q_simplices], Vec::new()), enumerate)
    .reduce(|| (vec![Vec::new(); number_q_simplices], Vec::new()), combine);
    for i in 0..number_q_simplices
    {
        result.0[i].push(i as Node);
    }   
    result
}

pub fn  get_q_graph_henri_improved<G: DirectedGraph + Sync + ?Sized + DirectedGraphNew+ AdjacencyMatrixGraph>(_graph: &G, _i:usize, _j:usize, _q:usize, number_q_simplices:usize, q_simplex_inclusions : &Vec<Vec<Node>>, simplicial_family: Vec<Vec<Node>>, _simplex_map: &HashMap<Simplex, Node>) -> EdgeListGraph {
    let enumerate = |mut q_graph: EdgeListGraph, q_simplex_ids: Vec<usize>| {
        for q_simplex_id in q_simplex_ids
        {
            let simplex_inclusion:&Vec<Node> = &q_simplex_inclusions[q_simplex_id as usize];
            
            let i_cobs = Vec::new();
            let j_cobs = Vec::new();
            for included_simplex in simplex_inclusion
            {
                    let _simplex:&Vec<Node> = &simplicial_family[*included_simplex as usize];
                    
                    //i_cobs.extend(cofaces::local_cofaces_id(graph, &simplex, simplex.len(), i as usize, &simplex_map));
                    //j_cobs.extend(cofaces::local_cofaces_id(graph, &simplex, simplex.len(), j as usize, &simplex_map));
            }
               
            for from in i_cobs
            {
                for to in &j_cobs
                {
                    if from != *to
                    {
                        q_graph.add_edge(from, *to);
                    }
                }
            }
        }
        q_graph
    };
        
    let combine = |mut right:  EdgeListGraph, mut left :  EdgeListGraph|
    {
        if left.nnodes > right.nnodes
        {
            for edge in right.edges()
            {
                left.add_edge(edge[0], edge[1]);
            }
            left
        }
        else
        {
            for edge in left.edges()
            {
                right.add_edge(edge[0], edge[1]);
            }
            right
        }
    };
    let q_graph = (0..number_q_simplices).into_par_iter().chunks(100).fold(|| EdgeListGraph::new_disconnected(0), enumerate).reduce(|| EdgeListGraph::new_disconnected(0), combine);
    q_graph
}

pub fn construct_q_graph_henri_improved(graph: EdgeMapGraph, flag_complex: Vec<Vec<Simplex>>, q:usize, i:usize, j:usize, simplex_map: &HashMap<Simplex, Node>, simplicial_family: Vec<Vec<Node>>) -> EdgeListGraph
{
    let (q_inclusions, inclusion_edges) = get_inclusion_with_q(&graph, i, j, q as usize, flag_complex[0].len(), &flag_complex, simplex_map);
    let mut q_graph = get_q_graph_henri_improved(&graph, i, j, q as usize,  flag_complex[0].len(), &q_inclusions, simplicial_family, simplex_map);
    for edge in inclusion_edges
    {
        q_graph.add_edge(edge[0], edge[1]);
    }
    dbg!(q_graph.edges().len());
    q_graph
}







// calculates the edges of the transitive reduction of the q-graph
fn construct_q_graph(graph: EdgeMapGraph, flag_complex: Vec<Vec<Simplex>>, q:Node, i:Node, j:Node) -> EdgeMapGraph
{
    let mut simplex_map:HashMap<Vec<Node>, usize> = HashMap::new();
    

    let f = |simplex: &Vec<Node>|{simplex.len() > q as usize};
    let simplicial_family:Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(f).collect();
    for (index, simplex) in simplicial_family.iter().enumerate()
    {
        simplex_map.insert(simplex.clone(), index);
    }
    let mut inclusion_graph = EdgeMapGraph::new_disconnected(simplicial_family.len());
    // inclusion edges
    let _max_dimension = flag_complex.len();
    for dimension in flag_complex.iter().rev()
    {
        if dimension[0].len() > (q+1) as usize
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
                            inclusion_graph.add_edge(*simplex_map.get(&face).unwrap() as Node, simplex_id as Node);
                        }
                    }
                }
            }
        }
        
    }
    // q-nearness edges of q+1-simplices
    let q_simplices = flag_complex.get((q) as usize).expect("No q simplex found");
    let mut q_near_graph = EdgeMapGraph::new_disconnected(simplicial_family.len());
    
    for simplex in q_simplices
    {
        let i_cofaces:Vec<Node> = crate::util::local_coboundaries(&graph, simplex,  i as usize);
        let j_cofaces:Vec<Node> = crate::util::local_coboundaries(&graph, simplex,  j as usize);

        for i_coboundary in i_cofaces
        {
            for j_coboundary in &j_cofaces
            {
                let mut simplex_i = simplex.clone();
                simplex_i.insert(i as usize, i_coboundary as Node);
                let mut simplex_j = simplex.clone();
                simplex_j.insert(j as usize, *j_coboundary as Node);
                let sigma = *simplex_map.get(&simplex_i).unwrap() as Node;
                let tau = *simplex_map.get(&simplex_j).unwrap() as Node;
                q_near_graph.add_edge(sigma, tau);
                // start a BFS to find all simplices that contain tau
                let tau_inclusion = bfs(&inclusion_graph, tau);
                let sigma_inclusion = bfs(&inclusion_graph, sigma);
                for from in sigma_inclusion
                {
                    for to in &tau_inclusion
                    {
                        if from != *to
                        {
                            q_near_graph.add_edge(from, *to);
                        }
                    }
                }
            }
        }
    }
    q_near_graph.add_graph_edges(inclusion_graph);
    q_near_graph
}

fn bfs(graph : &EdgeMapGraph, start : Node) -> Vec<Node>
{
    let mut visited = vec![false;graph.nnodes()];
    let mut queue = vec![start];
    let mut visited_nodes = vec![start];
    while queue.len() != 0 {
        let current_node = queue.pop().unwrap();
        visited[current_node as usize] = true;
        for node in graph.edges_from(current_node)
        {
            if !visited[node as usize]
            {
                visited_nodes.push(node);
                queue.push(node);
            }
            
        }
    }
    visited_nodes
}









