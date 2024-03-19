use crate::{graph::*,cofaces};
use libc::FALLOC_FL_UNSHARE_RANGE;
use rayon::prelude::*;
use std::sync::Mutex;
use std::{collections::HashMap, ops::Index, iter::Enumerate, hash::Hash, clone};
use crate::graph::{Node, Edge, EdgeMapGraph, DirectedGraphNew, AdjacencyMatrixGraph, EdgeListGraph};
use std::cmp::max;
use std::sync::atomic::AtomicU64;
use itertools::Itertools;
use itertools::enumerate;


pub fn compute_q_graph_fold_reduce(g:EdgeMapGraph, q:usize, i:usize, j:usize, flag_complex: &Vec<Vec<Simplex>>) -> EdgeListGraph
{
    let simplicial_family: Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(|simplex: &Vec<Node>|{simplex.len() > q as usize}).collect();
    //There is probably no way to parallelize this reasonably
    let mut simplex_map : HashMap<Simplex, Node> = HashMap::new();
    for (index, simplex) in enumerate(&simplicial_family)
    {
        simplex_map.insert(simplex.clone(), index as Node);
    }
    let (i_j_cofaces, (q_plus_one_inclusions, inclusion_edges))
    = rayon::join(
        ||{get_i_j_cofaces(&g, i, j, q, flag_complex[0].len(), &flag_complex[1][..], &simplex_map)},
        ||{get_inclusion(&g, i, j, q, flag_complex[1].len(), &flag_complex, &simplex_map)});
    //let simplex_count = flag_complex.iter().fold(0, |acc, x|{acc + x.len()});
    let mut q_graph = get_q_graph(&g, i, j, q, 0, flag_complex[0].len(), &i_j_cofaces, &q_plus_one_inclusions);
    dbg!(inclusion_edges.len());
    for edge in inclusion_edges
    {
        q_graph.add_edge(edge[0], edge[1]);
    }
    dbg!(q_graph.edges().len());
    q_graph
}





pub fn get_i_j_cofaces<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:usize, j:usize, q:usize, number_q_simplices:usize, q_plus_one_simplices: &[Vec<Node>], simplex_map: &HashMap<Simplex, Node>) -> Vec<(Vec<Node>, Vec<Node>)> {
    
    let mut enumerate = |mut result:  (Vec<(Vec<Node>, Vec<Node>)>, usize), simplex: &Vec<Node>| {

        let mut i_boundary :Vec<Node> = vec![0;simplex.len()];
            i_boundary =simplex.clone();
            i_boundary.remove(i as usize);
            let mut j_boundary :Vec<Node> = vec![0;simplex.len()];
            j_boundary = simplex.clone();
            j_boundary.remove(j as usize);
            
            let i_q_simplex = *simplex_map.get(&i_boundary).unwrap();
            let j_q_simplex = *simplex_map.get(&j_boundary).unwrap();
            let index = *simplex_map.get(simplex).unwrap();
            result.0[i_q_simplex as usize].0.push(index);
            result.0[j_q_simplex as usize].1.push(index);
            result.1 += 1;
        
        result

    };

    let combine = |mut right:  (Vec<(Vec<Node>, Vec<Node>)>, usize), mut left :  (Vec<(Vec<Node>, Vec<Node>)>, usize)|
    {
        let mut bigger;
        let smaller;
        if left.1 > right.1
        {
            bigger = left.0;
            smaller = right.0;
        }
        else
        {
            bigger = right.0;
            smaller = left.0;
        }
        let mut index = 0;
        for (i, j) in smaller
        {
            bigger[index].0.extend(i);
            bigger[index].1.extend(j);
            index += 1;
        }
        (bigger, left.1 + right.1)
    };


    q_plus_one_simplices.into_par_iter().fold(|| (vec![(Vec::new(), Vec::new()); number_q_simplices], 0), enumerate)
    .reduce(|| (vec![(Vec::new(), Vec::new()); number_q_simplices], 0), combine).0
}




pub fn get_inclusion<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:usize, j:usize, q:usize, number_q_plus_one_simplices:usize, flag_complex: &Vec<Vec<Simplex>>, simplex_map: &HashMap<Simplex, Node>) -> (Vec<Vec<Node>>, Vec<Edge>) {
    let mut enumerate = |mut result: (Vec<Vec<Node>>, Vec<Edge>), simplex: &Vec<Node>| {
        let simplex_id = *simplex_map.get(simplex).unwrap();
        for face_len in q as usize..simplex.len()
        {
            let faces : Vec<Vec<Node>> = simplex.clone().into_iter().combinations(face_len).collect();
            for face in faces
            {
                let face: Vec<Node> = face.into_iter().map(|x| x).collect();
                
                if face.len() > q as usize 
                {
                    let face_id = *simplex_map.get(&face).unwrap();
                    if face.len() == (q+2) as usize
                    {
                        result.0[(face_id - flag_complex[0].len() as Node) as usize].push(simplex_id);
                    }
                    result.1.push([face_id as Node, simplex_id as Node]);
                }
            }
        }
        (result.0, result.1)
    };
    
    let combine = |mut right:  (Vec<Vec<Node>>, Vec<Edge>), mut left :  (Vec<Vec<Node>>, Vec<Edge>)|
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
        bigger_edges.extend(smaller_edges);
        let number_cob_small = smaller_inclusion.iter().fold(0, |a, x| {a +x.len()});
        let number_cob_big = bigger_inclusion.iter().fold(0, |a, x| {a +x.len()});
        for index in 0..smaller_inclusion.len()
        {
            for x in smaller_inclusion[index].clone()
            {
                bigger_inclusion[index].push(x);
            }
        }
        let number_cob_after = bigger_inclusion.iter().fold(0, |a, x| {a +x.len()});

        assert!(number_cob_after == number_cob_small + number_cob_big);
        
        
        (bigger_inclusion, bigger_edges)
    };

    let mut result = flag_complex.into_par_iter().skip(1).flatten().fold(|| (vec![Vec::new(); number_q_plus_one_simplices], Vec::new()), enumerate)
    .reduce(|| (vec![Vec::new(); number_q_plus_one_simplices], Vec::new()), combine);
    for i in 0..number_q_plus_one_simplices
    {
        result.0[i].push(i as Node + flag_complex[0].len() as Node);
    }   
    result
}
    
pub fn get_q_graph<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:usize, j:usize, q:usize, number_simplices:usize, number_q_simplices:usize, i_j_cofaces: &Vec<(Vec<Node>, Vec<Node>)>, q_plus_one_simplex_inclusions : &Vec<Vec<Node>>) -> EdgeListGraph {
    let enumerate = |mut q_graph: EdgeListGraph, q_simplex_id: usize| {
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
                q_graph.add_edge(sigma, tau);    
                for from in sigma_inclusion
                {
                    for to in tau_inclusion
                    {
                        if *from != *to
                        {
                            q_graph.add_edge(*from, *to);
                        }
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
    let q_graph = (0..number_q_simplices).into_par_iter().fold(|| EdgeListGraph::new_disconnected(0), enumerate).reduce(|| EdgeListGraph::new_disconnected(0), combine);
    q_graph
}
