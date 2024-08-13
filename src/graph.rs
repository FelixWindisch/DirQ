pub type Node = u64;
pub type Edge = [Node; 2];
pub type Simplex = Vec<Node>;
use rand::prelude::*;
use serde::{Serialize, Deserialize};

use indexmap::set::IndexSet;
use std::cmp::{min, max};
use super::complex;




pub trait AdjacencyMatrixGraph
{
    //returns a row of the out_matrix as a slice
    fn edges_from_as_bitmap<'a>(&'a self, node:Node) -> &'a[Chunk];

    //returns a column of the out_matrix as a slice
    fn edges_to_as_bitmap<'a>(&'a self, node:Node) -> &'a[Chunk];

    fn row_len(&self) -> usize;

    fn print_adjacency_matrix(&self);

    fn transitive_closure(&self);

    fn add_graph_edges(&mut self, other:EdgeMapGraph);
}

impl AdjacencyMatrixGraph for EdgeMapGraph
{
    fn edges_from_as_bitmap<'a>(&'a self, node:Node) -> &'a[Chunk] {
        if node as usize >= self.nnodes {
            panic!("node out of bounds: {} >= {}", node, self.nnodes);
        }
        return &self.out_matrix[node as usize * self.row_len..node as usize *self.row_len + self.row_len];
    }

    fn edges_to_as_bitmap<'a>(&'a self, node:Node) -> &'a[Chunk] {
        if node as usize >= self.nnodes {
            panic!("node out of bounds: {} >= {}", node, self.nnodes);
        }
        return &self.out_matrix_transpose[node as usize * self.row_len..node as usize *self.row_len + self.row_len];
    }

    fn row_len(&self) -> usize{
        self.row_len
    }

    fn print_adjacency_matrix(&self)
    {
        for i in 0..self.nnodes
        {
            for j in 0..self.row_len
            {
                let x = self.out_matrix_transpose[i * self.row_len + j];
                let str_x = format!("{x:064b}");
                let reversed:Vec<u8> = str_x.as_bytes().iter().copied().rev().collect();
                for index in 0..64
                {
                    print!("{}", reversed[index]-48);
                }
            }
            println!("")
        }
    }

    //TODO implement as addition of binary matrices
    fn add_graph_edges(&mut self, other:EdgeMapGraph) {
        for edge in other.edges()
        {
            self.add_edge(edge[0], edge[1]);
        }
    }

    fn transitive_closure(&self) 
    {
        // TODO
    }
}

pub trait DirectedGraph: Sync {
    fn has_edge(&self, from: Node, to: Node) -> bool;
    fn add_edge(&mut self, from: Node, to: Node);
    fn remove_edge(&mut self, from: Node, to: Node);
    fn nnodes(&self) -> usize;
    fn iter_nodes(&self) -> std::ops::Range<Node> {
        (0.. (self.nnodes() as Node)).into_iter()
    }

    fn edges(&self) -> Vec<[Node; 2]> {
        let mut edges = vec![];
        for from in self.iter_nodes() {
            for to in self.iter_nodes() {
                if self.has_edge(from, to) {
                    edges.push([from, to]);
                }
            }
        }
        return edges
    }
    fn set_edge(&mut self, from: Node, to: Node, create: bool) {
        if create {
            self.add_edge(from, to);
        }
        else {
            self.remove_edge(from, to);
        }
    }

    

    fn edges_from(&self, from: Node) -> Vec<Node>{
        let mut res = vec![];
        for to in self.iter_nodes() {
            if self.has_edge(from, to) {
                res.push(to);
            }
        }
        return res;
    }

    fn edges_to(&self, to: Node) -> Vec<Node>{
        let mut res = vec![];
        for from in self.iter_nodes() {
            if self.has_edge(from, to) {
                res.push(from);
            }
        }
        return res;
    }

}

pub trait DirectedGraphNew: DirectedGraph + Sized {
    fn new_disconnected(nnodes: usize) -> Self;

    fn gen_seo_er<R: Rng>(nnodes: u64, p: f64, rng: &mut R) -> Self {
        let mut g = Self::new_disconnected(nnodes as usize);
        for i in 0..nnodes {
            for j in (i+1)..nnodes {
                if rng.gen::<f64>() < p {
                    if rng.gen::<f64>() < 0.5 {
                        g.add_edge(i,j);
                    } else {
                        g.add_edge(j,i);
                    }
                }
            }
        }
        return g;
    }

    fn subgraph<G: DirectedGraph>(ori: &G, vertices: &[Node]) -> Self {
        // macht komische dinge wenn Vertices in dem Slice doppelt vorkommen
        // FIXME?
        let mut sub = Self::new_disconnected(vertices.len());
        for (new_from, &ori_from) in (0..).zip(vertices.iter()) {
            for (new_to, &ori_to) in (0..).zip(vertices.iter()) {
                if ori.has_edge(ori_from, ori_to) {
                    sub.add_edge(new_from, new_to);
                }
            }
        }
        return sub;
    }
    //Should this be called clone()?
    fn copy<G: DirectedGraph>(g: &G) -> Self {
        let mut n = Self::new_disconnected(g.nnodes());
        for [a,b] in g.edges() {
            n.add_edge(a, b);
        }
        return n
    }
}


pub trait DirectedGraphExt: DirectedGraph {
    fn compute_maximal_cliques(&self) -> Vec<Vec<Node>> {
        // TODO: Kommentierung
        // undirected flagser statt eigenbau?
        let mut r = vec![];
        let x = vec![];
        let p = (0..(self.nnodes() as Node)).collect();
        let mut res = vec![];
        self.bron_kerbosch(&mut r, p, x, &mut res);
        return res;
    }
    fn compute_cliques(&self) -> Vec<Vec<Node>> {
        let mut edges = self.edges();
        for e in &mut edges {
            e.sort();
        }
        let mut normalized_graph = EdgeMapGraph::new_disconnected(self.nnodes());
        for [a,b] in edges {
            normalized_graph.add_edge(a, b);
        }
        let mut result = vec![];
        let mut add = |cell: &[Node]| {
            result.push(cell.into())
        };
        complex::for_each_cell(&normalized_graph, &mut add, 0, self.nnodes());
        return result
    }
    /// recursion for compute_maximal_cliques
    fn bron_kerbosch(&self, r: &mut Vec<Node>, p: Vec<Node>, mut x: Vec<Node>, res: &mut Vec<Vec<Node>>) {
        // from wikipedia: https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        // algorithm BronKerbosch1(R, P, X) is
        // if P and X are both empty then
        //     report R as a maximal clique
        // for each vertex v in P do
        //     BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        //     P := P \ {v}
        //     X := X ⋃ {v}
        
        if p.is_empty() && x.is_empty() {
            res.push(r.clone());
        }
        for (i,&v) in (0..).zip(&p) {
            let newp = p[i..].iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let newx = x.iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let addv = !r.contains(&v);
            if addv {
                r.push(v);
            }
            self.bron_kerbosch(r, newp, newx, res);
            if addv {
                r.pop();
            }
            x.push(v);
        }
    }

    fn undirected_edges(&self) -> Vec<Edge> {
        let mut undirected_edges = self.edges();
        for e in &mut undirected_edges {
            let a = max(e[0], e[1]);
            let b = min(e[0], e[1]);
            *e = [a, b];
        }
        undirected_edges.sort_unstable();
        undirected_edges.dedup();
        return undirected_edges;
    }

    fn flagser_count(&self) -> Vec<usize> {
        //crate::flagser::count_unweighted(self.nnodes(), &self.edges())
        complex::count_cells(self)
    }

}

impl<G: DirectedGraph> DirectedGraphExt for G {}

pub type Chunk = u64;
const CHUNK_SIZE: usize = Chunk::BITS as usize;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeMapGraph{
    nnodes: usize,
    edges: IndexSet<Edge>,
    double_edges: IndexSet<Edge>,
    //This is an adjacency matrix
    out_matrix: Vec<Chunk>,
    out_matrix_transpose: Vec<Chunk>,
    row_len: usize,
}

impl DirectedGraphNew for EdgeMapGraph {
    fn new_disconnected(nnodes: usize) -> Self {
        // nnodes / Chunk::BITS, rounded up.
        let row_len = (nnodes + Chunk::BITS as usize - 1) / (Chunk::BITS as usize);

        let matsize = row_len.checked_mul(nnodes).expect("size of adjacency matrix overflows");
        let out_matrix = vec![0; matsize];
        let out_matrix_transpose = vec![0; matsize];
        //let out_degrees = vec![0; nnodes];
        EdgeMapGraph 
        { 
            nnodes, 
            edges: Default::default(), 
            double_edges: Default::default(), 
            out_matrix, 
            out_matrix_transpose, 
            row_len
        }
    }
}

impl DirectedGraph for EdgeMapGraph {
    //TODO: Override default implementation for edges_to and edges_from to make use of out_matrix instead
    fn nnodes(&self) -> usize {
        self.nnodes
    }
    fn has_edge(&self, from: Node, to: Node) -> bool {
        if from as usize >= self.nnodes {
            panic!("from out of bounds: {} >= {}", from, self.nnodes);
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds: {} >= {}", to, self.nnodes);
        }
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        return self.out_matrix[row * self.row_len + toh] & (1<<tol) != 0;
    }
    fn add_edge(&mut self, from: Node, to: Node) {
        // update hashmaps
        if from != to
        {
            self.edges.insert([from, to]);
            if self.has_edge(to, from) {
                let big = max(from, to);
                let small = min(from, to);
                self.double_edges.insert([big, small]);
            }

            // update bits
            let toh = to as usize / CHUNK_SIZE;
            let tol = to as usize % CHUNK_SIZE;
            let row = from as usize;
            self.out_matrix[row * self.row_len + toh] |= 1<<tol;
            let toh_t = from as usize / CHUNK_SIZE;
            let tol_t = from as usize % CHUNK_SIZE;
            let row_t = to as usize;
            self.out_matrix_transpose[row_t * self.row_len + toh_t] |= 1<<tol_t;
        }
    }
    fn remove_edge(&mut self, from: Node, to: Node){
        // update hashmaps
        if self.has_edge(from, to) && self.has_edge(to, from) {
            let big = max(from, to);
            let small = min(from, to);
            self.double_edges.shift_remove(&[big, small]); // does nothing if it does not exist
        }
        self.edges.shift_remove(&[from, to]);

        //  update bits
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        self.out_matrix[row * self.row_len + toh] &= !(1<<tol);
        let toh_t = from as usize / CHUNK_SIZE;
        let tol_t = from as usize % CHUNK_SIZE;
        let row_t = to as usize;
        self.out_matrix_transpose[row_t * self.row_len + toh_t] &= !(1<<tol_t);
    }

    // FIXME: API ändern, dass nicht kopiert werden muss
    fn edges(&self) -> Vec<[Node; 2]> {
        self.edges.iter().cloned().collect()
    }
}

impl EdgeMapGraph {
    pub fn sample_edge<R: Rng>(&self, rng: &mut R) -> Option<Edge> {
        if self.edges.len() == 0 {
            return None;
        }
        let i = rng.gen_range(0 .. self.edges.len());
        return Some(self.edges[i]);
    }
    pub fn sample_double_edge<R: Rng>(&self, rng: &mut R) -> Option<Edge> {
        if self.double_edges.len() == 0 {
            return None;
        }
        let i = rng.gen_range(0 .. self.double_edges.len());
        return Some(self.double_edges[i]);
    }

}


/// Stores a graph as a list of edges. Space Complexity O(|E|).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeListGraph{
    pub nnodes: usize,
    edges: IndexSet<Edge>
}

impl DirectedGraph for EdgeListGraph 
{
    fn nnodes(&self) -> usize {
        self.nnodes
    }
    fn has_edge(&self, from: Node, to: Node) -> bool {
        for edge in &self.edges
        {
            if edge[0] == from && edge[1] == to
            {
                return true;
            }
        }
        return false;
    }
    fn add_edge(&mut self, from: Node, to: Node) {
        self.edges.insert([from, to]);
    }
    fn remove_edge(&mut self, from: Node, to: Node){
        self.edges.shift_remove(&[from, to]);
    }

    // FIXME: API ändern, dass nicht kopiert werden muss
    fn edges(&self) -> Vec<[Node; 2]> {
        self.edges.iter().cloned().collect()
    }
}


impl DirectedGraphNew for EdgeListGraph {
    fn new_disconnected(_nnodes: usize) -> Self 
    {
        EdgeListGraph 
        { 
            nnodes: 0, 
            edges: Default::default()
        }
    }
}