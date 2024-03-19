extern crate array_tool;
use std::collections::HashMap;

use bitvec::order::Msb0;

use bitvec::{vec::BitVec, *};
use crate::graph::Node;

pub use crate::graph::{Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = crate::graph::EdgeMapGraph;









/// DEPRECATED
/// Finds cofaces for a simplex tau as pairs of (index, vertex). 
/// The bit manipulations here are slow because I do not read them from the transition matrix.
pub fn local_cofaces_naive(g: &Graph, tau: &Vec<Node>, dimension: usize) -> Vec<(usize, usize)> {
    let mut co_b: Vec<(usize, usize)> = Vec::new();

    for i in 0..dimension {
        let mut c: BitVec<u8, Msb0> = bitvec![u8, Msb0; 1; 2000];
        for j in 0..i {
            let _node = tau[j];
            let mut out = bitvec![u8, Msb0; 0; 2000];
            for v in g.edges_from(tau[j]) {
                out.set(v as usize, true);
            }
            c = c & out;
        }
        for j in i..dimension - 1 {
            let _node = tau[j];
            let mut in_ = bitvec![u8, Msb0; 0; 2000];
            for v in g.edges_to(tau[j]) {
                in_.set(v as usize, true);
            }
            c = c & in_;
        }
        for v in 0..2000 {
            if *c.get(v).unwrap() {
                co_b.push((i, v))
            }
        }
    }
    co_b
}
