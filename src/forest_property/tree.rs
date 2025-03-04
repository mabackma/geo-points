use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Normal};
use wasm_bindgen::prelude::wasm_bindgen;
use web_sys::console::log_1;

#[derive(Default, Debug, Clone, Copy)]
#[wasm_bindgen]
pub struct Tree {
    stand_id: f64,
    species: u8,
    tree_height: f32,
    position: (f64, f64, f64),
    tree_status: f64,
}

impl Tree {
    pub fn new(
        stand_id: f64, 
        species: u8, 
        mean_height: f32, 
        position: (f64, f64, f64), 
        tree_status: Option<f64>
    ) -> Self {

        let height = calculate_height(mean_height);

        Tree {
            stand_id,
            species,
            tree_height: height,
            position,
            tree_status: tree_status.unwrap_or(1.0),
        }
    }

    pub fn species(&self) -> u8 {
        self.species
    }

    pub fn tree_height(&self) -> f32 {
        self.tree_height
    }

    pub fn position(&self) -> (f64, f64, f64) {
        self.position
    }

    pub fn tree_status(&self) -> f64 {
        self.tree_status
    }

    pub fn stand_id(&self) -> f64 {
        self.stand_id
    }

    pub fn cut_tree(&mut self) {
        self.tree_status = 0.0;
        log_1(&format!("Tree status after cutting: {}", self.tree_status).into());
    }

    pub fn set_position(&mut self, position: (f64, f64, f64)) {
        self.position = position;
    }
}

pub fn calculate_height(base_height: f32) -> f32 {
    
    // Create a Gaussian distribution with mean 1 and standard deviation 0.07
    let normal_dist = Normal::new(1.0, 0.07).unwrap();
    let mut rng = thread_rng();
    
    // Sample from the distribution to get the height ratio
    let mut height_ratio: f32 = normal_dist.sample(&mut rng);

    // Cap the height ratio between 0.9 and 1.1
    if height_ratio > 1.1 {
        height_ratio = rng.gen_range(0.9..=1.1);
    }

    base_height * height_ratio
}