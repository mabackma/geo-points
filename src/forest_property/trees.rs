use wasm_bindgen::prelude::wasm_bindgen;

use crate::forest_property::tree::Tree;

#[wasm_bindgen]
pub struct Trees {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    species: Vec<u8>,
    height: Vec<f32>,
    status: Vec<u8>,
    stand_number: Vec<u8>,
}

#[wasm_bindgen]
impl Trees {
    #[wasm_bindgen(constructor)]
    pub fn new(
        x: Vec<f32>,
        y: Vec<f32>,
        z: Vec<f32>,
        species: Vec<u8>,
        height: Vec<f32>,
        status: Vec<u8>,
        stand_number: Vec<u8>,
    ) -> Trees {
        Trees {
            x,
            y,
            z,
            species,
            height,
            status,
            stand_number,
        }
    }

    pub fn insert(&mut self, tree: Tree) {
        self.x.push(tree.position().0 as f32);
        self.y.push(tree.position().1 as f32);
        self.z.push(tree.position().2 as f32);
        self.species.push(tree.species());
        self.height.push(tree.tree_height());
        self.status.push(tree.tree_status() as u8);
        self.stand_number.push(tree.stand_number() as u8);
    }

    // Getters for each field
    #[wasm_bindgen(getter)]
    pub fn x(&self) -> Vec<f32> {
        self.x.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn y(&self) -> Vec<f32> {
        self.y.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn z(&self) -> Vec<f32> {
        self.z.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn specie(&self) -> Vec<u8> {
        self.species.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn height(&self) -> Vec<f32> {
        self.height.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn status(&self) -> Vec<u8> {
        self.status.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn stand_number(&self) -> Vec<u8> {
        self.stand_number.clone()
    }
}
