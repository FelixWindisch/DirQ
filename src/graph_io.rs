use std::fs::File;
use std::io::Read;

use std::collections::HashSet;
use crate::graph::*;


use super::graph;

pub fn read_flag_file<G: graph::DirectedGraphNew>(fname:&str) -> G {
    let mut file = File::open(fname).expect("could not find .flag input file");
    let mut fcontents = String::new();
    file.read_to_string(&mut fcontents).expect("File Parsing failed");

    let mut lines = fcontents.lines();
    lines.next(); // skip dim 0
    let nnodes = lines.next().expect("File Parsing failed").split(' ').filter(|s| *s != "").count(); // todo: mehr beleidigungen // todo: better parsing
    let mut graph = G::new_disconnected(nnodes);
    lines.next(); // skip dim 1
    for line in lines {
        let mut ijw = line.split(' ').filter(|s| *s != "");
        if let (Some(i), Some(j)) = (ijw.next(), ijw.next()) {
            graph.add_edge(i.parse().expect("File Parsing failed"), j.parse().expect("File Parsing failed"));
        }
    }
    return graph;
}

pub fn save_flag_file<G: graph::DirectedGraph>(fname:&str, graph:&G) {
    let mut nnodes:usize = 0;
    let mut edges = graph.edges();
    let mut edge_list = "".to_string();
    edges.sort_unstable();
    for [i,j] in edges {
        nnodes = std::cmp::max(nnodes, std::cmp::max( i as usize, j as usize));
        edge_list += &(format!("{} {} 1\n", i, j));
    }
    let mut content = "dim 0:\n".to_string();
    content += &("1 ".repeat(nnodes+1).to_owned() + "\n"); // add vertices
    content += "dim 1:\n";
    content += &edge_list;
    std::fs::write(fname, content).expect("Unable to write file");
}

pub fn save_flag_file_edge_list(fname:&str, edge_set: HashSet<(Simplex, Simplex)>) {
    let _nnodes:usize = 0;
    let mut edge_list = "".to_string();
    //edge_set.sort_unstable();
    for (sigma,tau) in edge_set {
        edge_list += &(format!("{:?} {:?}\n", sigma, tau));
    }
    std::fs::write(fname, edge_list).expect("Unable to write file");
}