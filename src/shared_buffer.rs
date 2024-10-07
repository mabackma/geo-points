use wasm_bindgen::prelude::*;
use web_sys::console::log_1;
use rand::Rng;

#[wasm_bindgen]
pub struct SharedBuffer {
    ptr: *mut f64,
    len: usize,
}

#[wasm_bindgen]
impl SharedBuffer {
    #[wasm_bindgen(constructor)]
    pub fn new(num_trees: usize) -> SharedBuffer {
        // Each tree will have 6 values: stand_number, x, y (f64), species (u8 stored as f64), tree_height (f32 stored as f64), tree_status (f64)
        let size = num_trees * 6;
        let buffer = vec![0f64; size].into_boxed_slice(); // Allocate memory
        let ptr = buffer.as_ptr() as *mut f64; // Get raw pointer to the buffer
        let len = buffer.len(); // Length of the buffer
        std::mem::forget(buffer); // Prevent Rust from freeing this memory
        SharedBuffer { ptr, len }
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
        if base + 5 < self.len / 6 && species != 0 {
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

    /// Marks trees by setting its status to 0.0 (stump).
    pub fn forest_clearing(&self, amount: usize, tree_count: usize) {
        let mut rng = rand::thread_rng();
        let mut index = 0;
        let mut indices = Vec::new();

        for _ in 0..amount {
            loop {
                index = rng.gen_range(0..tree_count);
                if !indices.contains(&index) {
                    indices.push(index);
                    break;
                }
            }

            let base = index * 6; 
            if base + 5 < self.len {
                unsafe {
                    *self.ptr.add(base + 5) = 0.0; // Change tree status to 0.0 (stump)
                }
            }
        }
    }

    // TODO: Implement forest thinning
    pub fn forest_thinning() {

    }
}


impl Drop for SharedBuffer {
    fn drop(&mut self) {
        unsafe {
            let _ = Vec::from_raw_parts(self.ptr, self.len, self.len);
        }
    }
}

