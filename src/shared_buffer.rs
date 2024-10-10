use wasm_bindgen::prelude::*;
use web_sys::console::log_1;
use rand::Rng;
use std::collections::HashSet;
use std::sync::Arc;

#[wasm_bindgen]
pub struct SharedBuffer {
    buffer: Arc<Vec<f64>>, // Use Arc to manage the buffer
    ptr: *mut f64,
    len: usize,
}

#[wasm_bindgen]
impl SharedBuffer {
    #[wasm_bindgen(constructor)]
    pub fn new(num_trees: usize) -> SharedBuffer {
        let size = num_trees * 6;
        let buffer = Arc::new(vec![0f64; size]); // Use Arc to manage buffer lifetime
        let ptr = buffer.as_ptr() as *mut f64; // Get raw pointer to the buffer
        let len = buffer.len();
        SharedBuffer { buffer, ptr, len }
    }

    pub fn ptr(&self) -> *mut f64 {
        self.ptr
    }

    pub fn set_ptr(&mut self, ptr: *mut f64) {
        self.ptr = ptr;
    }

    pub fn len(&self) -> usize {
        self.len
    }

    /// Fills the buffer with data for a single tree (stand_number, x, y, species, tree_height, tree_status=1.0)
    /// `index` is the index of the tree in the buffer (0-based)
    pub fn fill_tree(&self, index: usize, stand_number: f64, x: f64, y: f64, species: u8, tree_height: f32) {
        let base = index * 6; // 6 values per tree: stand_number, x, y, species, tree_height, tree_status
        if base + 6 < self.len && species != 0 {
            unsafe {
                *self.ptr.add(base) = stand_number;     // stand id in f64
                *self.ptr.add(base + 1) = x;           // x coordinate
                *self.ptr.add(base + 2) = y;       // y coordinate
                *self.ptr.add(base + 3) = species as f64; // species as u8 stored in f64
                *self.ptr.add(base + 4) = tree_height as f64;     // tree_height 
                *self.ptr.add(base + 5) = 1.0;     // tree_status (1.0 = tree, 0.0 = stump)
            }
        }
    }

    /// Clears trees from chosen stand. 
    pub fn forest_clearing(&self, stand_number: f64, amount: usize, tree_count: usize, area_ratio: f64) {
        let mut rng = rand::thread_rng();
        let mut indices = HashSet::new();
        let mut trees_cleared = 0;
        let trees_to_cut = (amount as f64 * area_ratio).floor() as usize;
    
        while trees_cleared < trees_to_cut {
            let index = rng.gen_range(0..tree_count);
            
            if indices.insert(index) { // Only process if index wasn't already processed
                let base = index * 6;

                unsafe {
                    if base + 6 < self.len && *self.ptr.add(base) == stand_number {
                        *self.ptr.add(base + 5) = 0.0; // Change tree status to 0.0 (stump)
                        trees_cleared += 1;
                    }
                }
            }
        }
    }

    // TODO: Implement forest thinning
    pub fn forest_thinning() {

    }

    pub fn log_buffer(&self) {
        unsafe {
            for i in 0..self.len {
                log_1(&format!("Buffer[{}]: {}", i, *self.ptr.add(i)).into());
            }
        }
    }
}

