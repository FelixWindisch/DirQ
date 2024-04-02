use crate::graph::*;

use rayon::prelude::*;
use std::collections::HashMap;
use crate::graph::{Node, Edge};
use std::cmp::max;

use itertools::enumerate;
// This is a translation of parts from the "flagser" C++ code of Daniel Lütgehetmann.
// Published as "Computing Persistent Homology of Directed Flag Complexes" https://doi.org/10.3390/a13010019
// https://github.com/luetge/flagser/blob/master/include/complex/directed_flag_complex.h


// Returns the directed flag complex of G as a Vec<Vec<Vec<Node>>>,
// as well as a Hashmap of simplices to indices and i/j-cofaces for all q-simplices
pub fn get_directed_flag_complex<G: DirectedGraph + Sync + ?Sized>(graph: &G, _i:usize, _j:usize, q:usize) -> Vec<Vec<Vec<Node>>> {
    
    let enumerate = |result: &mut Vec<Vec<Simplex>>, simplex: &[Node]| {
        // Should this closure be declared outside?
        //if x.len() > q{
            let mut x: Vec<Node> = vec![];
            x.extend_from_slice(simplex);
            if x.len() > q
            {
                if result.len() < simplex.len()-q {
                    result.push(vec![x]);
                } else {
                    result[simplex.len() - 1 - q].push(x);
                }
            }
    };

    let combine = |left:  Vec<Vec<Simplex>>, right: Vec<Vec<Simplex>>| 
    {
        let mut flag_complex_l = left;
        let mut flag_complex_r = right;
        
        let max_dimension = max(flag_complex_l.len(), flag_complex_r.len());
        while max_dimension > flag_complex_l.len()
        {
            flag_complex_l.push(Vec::new());
        }
        while max_dimension > flag_complex_r.len()
        {
            flag_complex_r.push(Vec::new());
        }
        for dimension in 0..max_dimension
        {
            flag_complex_l[dimension].extend(flag_complex_r[dimension].clone());
        }        
        
        flag_complex_l
    };

    crate::complex::for_each_cell_par(graph, &enumerate, 0, graph.nnodes()).reduce(||{Vec::new()}, &combine)
}

pub fn get_simplex_map(flag_complex: &Vec<Vec<Simplex>>, q:usize) -> HashMap<Simplex, Node>
{
    //.filter(|simplex: &Vec<Node>|{simplex.len() > q as usize})
    let simplicial_family: Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().collect();
    let mut simplex_map : HashMap<Simplex, Node> = HashMap::new();
    for (index, simplex) in enumerate(&simplicial_family)
    {
        simplex_map.insert(simplex.clone(), index as Node);
    }
    simplex_map
}



// Returns the directed flag complex of G as a Vec<Vec<Vec<Node>>>,
// as well as a Hashmap of simplices to indices and i/j-cofaces for all q-simplices
pub fn enumerate_cells_with_q_nearness<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:Node, j:Node, q:usize) -> (Vec<Vec<Vec<Node>>>, HashMap<Simplex, Node>, Vec<(Vec<Node>, Vec<Node>)>) {
    let mut flag_complex: Vec<Vec<Vec<Node>>> = Vec::new();
    let mut i_j_cofaces: Vec<(Vec<Node>, Vec<Node>)> = Vec::new();
    let mut simplex_map: HashMap<Simplex, Node> = HashMap::new();
    let mut n_simplices = 0;
    
    let mut enumerate = |simplex: &[Node]| {
        // Should this closure be declared outside?
        let mut get_index_or_add_simplex = |x: Vec<Node>|
        {
            if simplex_map.contains_key(&x)
            {
                *simplex_map.get(&x).unwrap()
            }
            else
            {
                simplex_map.insert(x, n_simplices);
                n_simplices += 1; 
                i_j_cofaces.push((Vec::new(), Vec::new())); 
                n_simplices - 1
            }
        };

        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(simplex);
        let mut index=0;
        if simplex.len() > q
        {
            index = get_index_or_add_simplex(x.clone());
        }
        // find the q+1-simplices
        if simplex.len() == q+2
        {
            let mut i_coboundary :Vec<Node> = vec![0;simplex.len()];
            i_coboundary.copy_from_slice(simplex);
            i_coboundary.remove(i as usize);
            let mut j_coboundary :Vec<Node> = vec![0;simplex.len()];
            j_coboundary.copy_from_slice(simplex);
            j_coboundary.remove(j as usize);
            
            let i_q_simplex = get_index_or_add_simplex(i_coboundary);
            let j_q_simplex = get_index_or_add_simplex(j_coboundary);
            i_j_cofaces[i_q_simplex as usize].0.push(index);
            i_j_cofaces[j_q_simplex as usize].1.push(index);
        }
        
        if flag_complex.len() < simplex.len() {
            flag_complex.push(vec![x]);
        } else {
            flag_complex[simplex.len() - 1].push(x);
        }
    };
    for_each_cell(graph, &mut enumerate, 0, graph.nnodes());
    return (flag_complex, simplex_map, i_j_cofaces);
}

// Returns the directed flag complex of G as a Vec<Vec<Vec<Node>>>,
// as well as a Hashmap of simplices to indices and i/j-cofaces for all q-simplices
pub fn enumerate_cells_with_q_nearness_par<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:Node, j:Node, q:usize) -> (Vec<Vec<Vec<Node>>>, HashMap<Simplex, Node>, Vec<(Vec<Node>, Vec<Node>)>) {
    
    let enumerate = |result: &mut( Vec<Vec<Simplex>>, HashMap<Simplex, Node>, Vec<(Vec<Node>, Vec<Node>)>), simplex: &[Node]| {
        // Should this closure be declared outside?
        let mut get_index_or_add_simplex = |x: Vec<Node>|
        {
            if result.1.contains_key(&x)
            {
                *result.1.get(&x).unwrap()
            }
            else
            {
                result.1.insert(x, result.1.len() as Node);
                result.2.push((Vec::new(), Vec::new())); 
                (result.1.len() - 1) as Node
            }
        };

        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(simplex);
        let mut index=0;
        if simplex.len() > q
        {
            index = get_index_or_add_simplex(x.clone());
        }
        // find the q+1-simplices
        if simplex.len() == q+2
        {
            let mut i_coboundary :Vec<Node> = vec![0;simplex.len()];
            i_coboundary.copy_from_slice(simplex);
            i_coboundary.remove(i as usize);
            let mut j_coboundary :Vec<Node> = vec![0;simplex.len()];
            j_coboundary.copy_from_slice(simplex);
            j_coboundary.remove(j as usize);
            
            let i_q_simplex = get_index_or_add_simplex(i_coboundary);
            let j_q_simplex = get_index_or_add_simplex(j_coboundary);
            result.2[i_q_simplex as usize].0.push(index);
            result.2[j_q_simplex as usize].1.push(index);
        }
        //if x.len() > q{
        if result.0.len() < simplex.len() {
            result.0.push(vec![x]);
        } else {
            result.0[simplex.len() - 1].push(x);
        }
        //}
    };

    let combine = |left: ( Vec<Vec<Simplex>>, HashMap<Simplex, Node>, Vec<(Vec<Node>, Vec<Node>)>), right: ( Vec<Vec<Simplex>>, HashMap<Simplex, Node>, Vec<(Vec<Node>, Vec<Node>)>)| 
    {
        let (mut flag_complex_l, mut simplex_map_l, i_j_cofaces_l) = left;
        let (mut flag_complex_r, simplex_map_r, _i_j_cofaces_r) = right;
        
        let max_dimension = max(flag_complex_l.len(), flag_complex_r.len());
        while max_dimension > flag_complex_l.len()
        {
            flag_complex_l.push(Vec::new());
        }
        while max_dimension > flag_complex_r.len()
        {
            flag_complex_r.push(Vec::new());
        }
        for dimension in 0..max_dimension
        {
            flag_complex_l[dimension].extend(flag_complex_r[dimension].clone());
        }
        simplex_map_l.extend(simplex_map_r);
        
        


        (flag_complex_l, simplex_map_l, i_j_cofaces_l)
    };

    for_each_cell_par(graph, &enumerate, 0, graph.nnodes()).reduce(||{(Vec::new(), HashMap::new(), Vec::new())}, &combine)
    //(Vec::new(), HashMap::new(), Vec::new())

}









pub struct FullEnumerationResult
{
    pub size:Node,
    pub flag_complex  : Vec<Vec<Vec<Node>>>,
    pub simplex_map: HashMap<Simplex, Node>,
    pub i_j_cofaces: Vec<(Vec<Node>, Vec<Node>)>,
    pub q_plus_one_simplex_inclusions:Vec<Vec<Node>>,
    pub inclusion_edges: Vec<Edge>
}

impl Default for FullEnumerationResult
{
    fn default() -> Self {
        FullEnumerationResult
        {
            size: 0,
            flag_complex: Vec::new(),
            simplex_map: HashMap::new(),
            i_j_cofaces: Vec::new(),
            q_plus_one_simplex_inclusions: Vec::new(),
            inclusion_edges: Vec::new()
        }
    }
}

impl FullEnumerationResult
{
    fn combine(&mut self, b: &mut FullEnumerationResult) -> &mut FullEnumerationResult
    {
        let smaller:  &mut FullEnumerationResult;
        let bigger:  &mut FullEnumerationResult;
        if self.size > b.size
        {
            smaller = b;
            bigger = self;
        }
        else
        {
            smaller = b;
            bigger = self;
        }
        
        bigger.size = bigger.size + smaller.size;
        //combine flag complexes
        let max_dimension = max(bigger.flag_complex.len(), smaller.flag_complex.len());
        while max_dimension > bigger.flag_complex.len()
        {
            bigger.flag_complex.push(Vec::new());
        }
        while max_dimension > smaller.flag_complex.len()
        {
            smaller.flag_complex.push(Vec::new());
        }
        for dimension in 0..max_dimension
        {
            bigger.flag_complex[dimension].extend(smaller.flag_complex[dimension].clone());
        }
        for entry in smaller.simplex_map.iter()
        {
            bigger.simplex_map.insert(entry.0.to_owned(), entry.1 + bigger.size);
        }
        //bigger.simplex_map.extend(smaller.simplex_map.clone());

        //This is really bad for performance
        for entry in smaller.i_j_cofaces.clone()
        {
            bigger.i_j_cofaces.push((entry.0.iter().map(|x| x+bigger.size).collect(), entry.0.iter().map(|x| x+bigger.size).collect()) );
        }
        bigger.i_j_cofaces.extend(smaller.i_j_cofaces.clone());
        

        for (index, entry) in smaller.q_plus_one_simplex_inclusions.iter().enumerate()
        {
            if entry.len() > 0
            {
                bigger.q_plus_one_simplex_inclusions[index].extend(entry.clone());
            }
        }

        bigger.inclusion_edges.extend(smaller.inclusion_edges.clone());

        bigger
    }
}


// Returns the directed flag complex of G as a Vec<Vec<Vec<Node>>>,
// as well as a Hashmap of simplices to indices and i/j-cofaces for all q-simplices
pub fn enumerate_cells_with_q_nearness_par_full<G: DirectedGraph + Sync + ?Sized>(graph: &G, i:Node, j:Node, q:usize) -> FullEnumerationResult {
    
    let enumerate = |result: &mut FullEnumerationResult, simplex: &[Node]| {
        // Should this closure be declared outside?
        let mut get_index_or_add_simplex = |x: Vec<Node>|
        {
            if result.simplex_map.contains_key(&x)
            {
                *result.simplex_map.get(&x).unwrap()
            }
            else
            {
                result.size += 1;
                result.simplex_map.insert(x, result.simplex_map.len() as Node);
                result.i_j_cofaces.push((Vec::new(), Vec::new())); 
                (result.simplex_map.len() - 1) as Node
            }
        };

        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(simplex);
        let mut index=0;
        if simplex.len() > q
        {
            index = get_index_or_add_simplex(x.clone());
        }
        // find the q+1-simplices
        if simplex.len() == q+2
        {
            let mut i_boundary :Vec<Node> = vec![0;simplex.len()];
            i_boundary.copy_from_slice(simplex);
            i_boundary.remove(i as usize);
            let mut j_boundary :Vec<Node> = vec![0;simplex.len()];
            j_boundary.copy_from_slice(simplex);
            j_boundary.remove(j as usize);
            
            let i_q_simplex = get_index_or_add_simplex(i_boundary);
            let j_q_simplex = get_index_or_add_simplex(j_boundary);
            result.i_j_cofaces[i_q_simplex as usize].0.push(index);
            result.i_j_cofaces[j_q_simplex as usize].1.push(index);
        }
        //if x.len() > q{
        if result.flag_complex.len() < simplex.len() {
            result.flag_complex.push(vec![x]);
        } else {
            result.flag_complex[simplex.len() - 1].push(x);
        }

        //}
    };

    let combine = |mut left: FullEnumerationResult, mut right: FullEnumerationResult| 
    {
        left.combine(&mut right);
        left
    };

    for_each_cell_par(graph, &enumerate, 0, graph.nnodes()).reduce(||{FullEnumerationResult::default()}, &combine)
    //(Vec::new(), HashMap::new(), Vec::new())

}








// Generates a Directed Flag Complex. result[2][1][0] contains the first node of the second simplex of dimension 2.
pub fn enumerate_cells<G: DirectedGraph + Sync + ?Sized>(graph: &G, q: usize) -> Vec<Vec<Vec<Node>>> {
    let mut result: Vec<Vec<Vec<Node>>> = Vec::new();
    let mut enumerate = |cell: &[Node]| {
        let mut x: Vec<Node> = vec![];
        x.extend_from_slice(cell);
        if x.len() > q
        {
            if result.len() < cell.len()-q {
                result.push(vec![x]);
            } else {
                result[cell.len() - 1 - q].push(x);
            }
        }
    };
    for_each_cell(graph, &mut enumerate, 0, graph.nnodes());
    return result;
}

pub fn for_each_cell_par<'a, G: DirectedGraph+Sync+?Sized, R, F>
    (graph: &'a G, f: F, min_dimension: usize, max_dimension: usize) -> impl ParallelIterator<Item=R> + 'a
where
    F: Fn(&mut R, &[Node]) + Send + Sync + 'a,
    R: Default+Send+Sync + 'a,
{
    (0 .. graph.nnodes() as Node).into_par_iter().fold(R::default, move |mut state, v| {
        let possible_next_vertices = graph.edges_from(v);
        let mut prefix = vec![v];
        do_for_each_cell(graph, &mut |simplex| f(&mut state, simplex), min_dimension, max_dimension, &mut prefix, &possible_next_vertices);
        state
    })
}

pub fn for_each_cell<G,F>
    (graph: &G, f: &mut F, min_dimension: usize, max_dimension: usize)
where
    G: DirectedGraph+?Sized,
    F: FnMut(&[Node]),
{
    let mut prefix = vec![];
    for v in graph.iter_nodes() {
        prefix.resize(1, v);
        prefix[0] = v;
        let possible_next_vertices = graph.edges_from(v);
        do_for_each_cell(graph, f, min_dimension, max_dimension, &mut prefix, &possible_next_vertices);
    }
}

pub fn for_each_cell_with_inclusion<G,F>
    (graph: &G, f: &mut F, min_dimension: usize, max_dimension: usize)
where
    G: DirectedGraph+?Sized,
    F: FnMut(&[Node], Node),
{
    let mut prefix = vec![];
    for v in graph.iter_nodes() {
        prefix.resize(1, v);
        prefix[0] = v;
        let possible_next_vertices = graph.edges_from(v);
        do_for_each_cell_with_inclusion(graph, f, min_dimension, max_dimension, &mut prefix, &possible_next_vertices);
    }
}

fn do_for_each_cell_with_inclusion<F,G>(graph: &G, f: &mut F, min_dimension: usize,
    max_dimension: usize,
    prefix: &mut Vec<Node>,
    possible_next_vertices: &[Node]
    )
where G: DirectedGraph+?Sized, F: FnMut(&[Node], Node)
{

    // Why wrapping_add???
    if prefix.len() == max_dimension.wrapping_add(1) {
        return;
    }

    for &vertex in possible_next_vertices {
        let new_possible_next_vertices = if prefix.len() == 0 {
            graph.edges_from(vertex)
        } else {
            possible_next_vertices.iter().cloned().filter(|&v| graph.has_edge(vertex, v)).collect()
        };
        
        prefix.push(vertex);
        if prefix.len() >= min_dimension + 1 {
            f(&prefix, vertex);
        }
            do_for_each_cell_with_inclusion(graph, f, min_dimension, max_dimension, prefix, &new_possible_next_vertices);
        prefix.pop();
}
}



fn do_for_each_cell<F,G>(graph: &G, f: &mut F, min_dimension: usize,
                         max_dimension: usize,
                         prefix: &mut Vec<Node>,
                         possible_next_vertices: &[Node]
                         )
    where G: DirectedGraph+?Sized, F: FnMut(&[Node])
{
    if prefix.len() >= min_dimension + 1 {
        f(&prefix);
    }
    // Why wrapping_add???
    if prefix.len() == max_dimension.wrapping_add(1) {
        return;
    }

    for &vertex in possible_next_vertices {
        let new_possible_next_vertices = if prefix.len() == 0 {
            graph.edges_from(vertex)
        } else {
            possible_next_vertices.iter().cloned().filter(|&v| graph.has_edge(vertex, v)).collect()
        };
        prefix.push(vertex);
        do_for_each_cell(graph, f, min_dimension, max_dimension, prefix, &new_possible_next_vertices);
        prefix.pop();
    }
}

pub fn count_cells_par<G: DirectedGraph+Sync+?Sized>(graph: &G) -> Vec<usize>{
    let count = |result: &mut Vec<usize>, cell: &[Node]| {
        if result.len() < cell.len() {
            result.resize(cell.len(), 0);
        }
        result[cell.len()-1] += 1;
    };
    let combine = |mut left: Vec<usize>, mut right: Vec<usize>| {
        if left.len() < right.len() {
            std::mem::swap(&mut left, &mut right)
        }
        for i in 0..right.len() {
            left[i] += right[i];
        }
        left
    };
    let res = for_each_cell_par(graph, &count, 0, graph.nnodes()).reduce(Vec::new, &combine);
    return res;
}

pub fn count_cells<G: DirectedGraph+Sync+?Sized>(graph: &G) -> Vec<usize>{
    let mut result = vec![];
    let mut count = |cell: &[Node]| {
        if result.len() < cell.len() {
            result.resize(cell.len(), 0);
        }
        result[cell.len()-1] += 1;
    };
    for_each_cell(graph, &mut count, 0, graph.nnodes());
    return result;
}


