

use eframe::{egui};
use tinyfiledialogs::{open_file_dialog,select_folder_dialog};

use crate::graph;

use crate::directed_q;
pub use graph::{EdgeMapGraph, EdgeListGraph, DirectedGraph, DirectedGraphExt, DirectedGraphNew, AdjacencyMatrixGraph};
pub type Graph = graph::EdgeMapGraph;





#[derive(PartialEq)]
#[derive(Eq)]
enum Definition { New, Original }

impl std::default::Default for Definition
{
    fn default() -> Self {
        return Definition::New
    }
}

impl std::fmt::Display for Definition
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self
        {
            Definition::New => write!(f, "new"),
            Definition::Original => write!(f, "original")
        }
        
    }
}

#[derive(PartialEq)]
enum Parallelization { SingleThread, Mutex, SplitAndMerge, BottomUp, Naive }

impl std::default::Default for Parallelization
{
    fn default() -> Self {
        return Parallelization::SingleThread
    }
}

#[derive(Default)]
pub struct DirQ 
{
    
    file : String, 
    out_folder : String, 
    q : String,
    q_val : u32,
    i : String,
    j : String,
    maximal_dimension : String,
    def : Definition,
    par: Parallelization,
    running : bool,
    finished : bool,
    written_file: String
}

impl DirQ {
    fn name() -> &'static str {
        "Fast Directed Q-Analysis"
    }
}


impl eframe::App for DirQ {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_pixels_per_point(1.5);

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.set_enabled(!self.running);

            ui.heading("DirQ: Fast Directed Q-Analysis");
            if ui.button("Input File").clicked() {
                self.file = open_file_dialog("Select input .flag file", "", None).unwrap_or(String::from(""));
            };
            ui.label(&self.file);
            ui.horizontal(|ui| {
                ui.label("q");
                ui.add_sized([6.,20.], egui::TextEdit::singleline(&mut self.q));
                ui.label("i");
                ui.add_sized([6.,20.], egui::TextEdit::singleline(&mut self.i));
                ui.label("j");
                ui.add_sized([6.,20.], egui::TextEdit::singleline(&mut self.j));
                
            
            });
            ui.horizontal(|ui| {
                ui.label("Maximal Simplex Dimension");
                ui.add_sized([6.,20.], egui::TextEdit::singleline(&mut self.maximal_dimension));
            });
            self.q_val = self.q.parse().unwrap_or_else(|_x| {self.q = String::from(""); 0});
            if self.q_val > 10
            {
                ui.label("The value for q is very high");
            }
            
            if self.i.parse().unwrap_or_else(|_x| {self.i = String::from(""); 0}) > self.q_val+1 && matches!(self.def, Definition::New)
            {
                ui.label("The value of i should not exceed q+1");
            }
            if self.j.parse().unwrap_or_else(|_x| {self.j = String::from(""); 0}) > self.q_val+1 && matches!(self.def, Definition::New)
            {
                ui.label("The value of j should not exceed q+1");
            }
            self.maximal_dimension.parse().unwrap_or_else(|_x| {self.maximal_dimension = String::from(""); 0}); 
            
            ui.horizontal(|ui| {
                ui.label("Definition:");
                ui.selectable_value(&mut self.def, Definition::New, "New");
                ui.selectable_value(&mut self.def, Definition::Original, "Original");
            });

            ui.horizontal(|ui| {
                ui.label("Parallelization:");
                ui.selectable_value(&mut self.par, Parallelization::SingleThread, "Single Threaded");
                ui.selectable_value(&mut self.par, Parallelization::Mutex, "Mutex");
                ui.selectable_value(&mut self.par, Parallelization::SplitAndMerge, "Split and Merge");
                ui.selectable_value(&mut self.par, Parallelization::BottomUp, "Bottom Up");
                ui.selectable_value(&mut self.par, Parallelization::Naive, "Naive");
            });
            if matches!(self.par, Parallelization::Naive)
            {
                ui.label("The naive algorithm is significantly slower");
            }
            if matches!(self.par, Parallelization::BottomUp)
            {
                ui.label("The bottom up algorithm currently does not support flagser file export and will export to a custom format");
            }
            
            if ui.button("Output Folder").clicked() {
                self.out_folder = select_folder_dialog("Select output folder", "").unwrap_or(String::from(""));
            };
            ui.label(&self.out_folder);
            ui.set_enabled(!self.q.is_empty() && !self.i.is_empty() && !self.j.is_empty()&& !self.file.is_empty() && !self.running);
            
            if ui.button("Run").clicked() {
                self.finished = false;
                self.running = true;
                let args = directed_q::Args
                {
                    max_dimension : self.maximal_dimension.parse().unwrap_or(100000000),
                    i : self.i.parse().unwrap_or_default(),
                    j : self.j.parse().unwrap_or_default(),
                    q : self.q.parse().unwrap_or_default(),
                input_file : self.file.clone(),
                output_folder : self.out_folder.clone(),
                parallelization : String::from(match &self.par {
                    Parallelization::SingleThread => "none", 
                    Parallelization::Mutex => "mutex",
                    Parallelization::SplitAndMerge => "split_and_merge", 
                    Parallelization::BottomUp => "bottom_up",
                    Parallelization::Naive => "naive",
                }),
                definition : String::from(match &self.def
                {
                    Definition::New => "new",
                    Definition::Original => "original", 
                }),
                benchmark : "".to_string(),
            };
            
                directed_q::run(args);
                self.running = false;
                self.finished = true;
            };
            ui.set_enabled(true);
            if ui.button("Quit").clicked() {
                std::process::exit(0);
            };
            self.written_file = "".to_string();
            self.written_file.push_str(&self.file.split('.').next().unwrap());
            self.written_file.push_str( &(format!("_q{}i{}j{}_{}", self.q, self.i, self.j, self.def)));
            if matches!(self.par, Parallelization::BottomUp)
            {
                self.written_file.push_str(".q_graph");
            }
            else
            {
                self.written_file.push_str(".flag");
            }
            if self.finished
            {
                ui.label("output written to ".to_string() + &self.written_file);
            }
        });
    }
}

pub fn do_gui() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions {
        ..eframe::NativeOptions::default()
    };

    eframe::run_native(
        DirQ::name(),
        native_options,
        Box::new(|_| Box::<DirQ>::default()),
    )
} 


