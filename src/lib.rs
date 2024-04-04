use graph::{EdgeMapGraph, Node};
use pyo3::prelude::*;


use pyo3::types::PyList;
use std::collections::HashMap;
use itertools::enumerate;
use numpy::PyReadonlyArray2;
use pyo3::{pymodule, types::PyModule, PyResult, Python};
pub mod graph;
pub mod complex;
pub mod graph_io;
pub mod directed_q;
pub mod util;
use directed_q::new::split_and_merge;
pub use graph::{Simplex, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};


const CHUNKSIZE:usize=1000000;

#[pyfunction]
fn compute_q_near_graph_old<'py>(py: Python<'py>, a: PyReadonlyArray2<usize>, q:usize, i:usize, j:usize) -> PyResult<&'py PyList>
{
    assert_eq!(a.shape()[0], a.shape()[1]);
    let graph_size = a.shape()[0];
    let mut g =  graph::EdgeMapGraph::new_disconnected(graph_size);
    for x in 0..graph_size
    {
        for y in 0..graph_size
        {
            if *a.get([x,y]).unwrap() == 1
            {   
                assert_ne!(x, y);
                 g.add_edge(x as Node, y as Node);
            }
        }   
    }
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplicial_family: Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(|simplex: &Vec<Node>|{simplex.len() > q as usize}).collect();
    let mut simplex_map : HashMap<Simplex, Node> = HashMap::new();
    for (index, simplex) in enumerate(&simplicial_family)
    {
        simplex_map.insert(simplex.clone(), index as Node);
    }
    let r = directed_q::old::single_thread::get_q_digraph(q, i, j, &flag_complex,  &simplex_map);

    let result = PyList::new(py, r.edges());
    Ok(result)
}

#[pyfunction]
fn compute_q_near_graph_new<'py>(py: Python<'py>, a: PyReadonlyArray2<usize>, q:usize, i:usize, j:usize ) -> PyResult<&'py PyList>
{
    assert_eq!(a.shape()[0], a.shape()[1]);
    let graph_size = a.shape()[0];
    let mut g =  graph::EdgeMapGraph::new_disconnected(graph_size);
    for x in 0..graph_size
    {
        for y in 0..graph_size
        {
            if *a.get([x,y]).unwrap() == 1
            {   
                assert_ne!(x, y);
                 g.add_edge(x as Node, y as Node);
            }
        }   
    }
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);
    let r = split_and_merge::get_q_digraph( q, i, j, &flag_complex, &simplex_map).edges();
    let result = PyList::new(py, r);
    Ok(result)
}

#[pyfunction]
fn compute_q_near_graph_old_from_file<'py>(py: Python<'py>, input_file:String, q:usize, i:usize, j:usize) -> PyResult<&'py PyList>
{
    let g : EdgeMapGraph = graph_io::read_flag_file(&input_file);
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplicial_family: Vec<Vec<Node>> = flag_complex.clone().into_iter().flatten().filter(|simplex: &Vec<Node>|{simplex.len() > q as usize}).collect();
    let mut simplex_map : HashMap<Simplex, Node> = HashMap::new();
    for (index, simplex) in enumerate(&simplicial_family)
    {
        simplex_map.insert(simplex.clone(), index as Node);
    }
    let r = directed_q::old::single_thread::get_q_digraph(q, i, j, &flag_complex, &simplex_map);

    let result = PyList::new(py, r.edges());
    Ok(result)
}

#[pyfunction]
fn compute_q_near_graph_new_from_file<'py>(py: Python<'py>, input_file:String, i:usize, j:usize, q:usize) -> PyResult<&'py PyList>
{
    let g : EdgeMapGraph = graph_io::read_flag_file(&input_file);
    let flag_complex = complex::get_directed_flag_complex(&g, i, j, q);
    let simplex_map = complex::get_simplex_map(&flag_complex, q);
    let r = split_and_merge::get_q_digraph( q, i, j, &flag_complex, &simplex_map).edges();
    let result = PyList::new(py, r);
    Ok(result)
}


#[pymodule]
fn dir_q_lib(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_q_near_graph_old, m)?)?;
    m.add_function(wrap_pyfunction!(compute_q_near_graph_new, m)?)?;
    m.add_function(wrap_pyfunction!(compute_q_near_graph_old_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(compute_q_near_graph_new_from_file, m)?)?;                                                                                                                                                                                                                                                                                                                                    
    Ok(())
}
