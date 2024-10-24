use crate::forest_property::tree::Tree;

use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::JsValue;
use web_sys::js_sys::{self, Array};
use js_sys::Function;

#[wasm_bindgen]
pub struct Trees {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    species: Vec<u8>,
    height: Vec<f32>,
    status: Vec<u8>,
    stand_number: Vec<u16>,
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
        stand_number: Vec<u16>,
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
        self.stand_number.push(tree.stand_number() as u16);
    }

    pub fn for_each(&self, callback: &Function) {
        for i in 0..self.x.len() {
            // Create a JS Array for the arguments
            let tree_args = Array::new();
            tree_args.push(&JsValue::from_f64(self.x[i] as f64));
            tree_args.push(&JsValue::from_f64(self.y[i] as f64));
            tree_args.push(&JsValue::from_f64(self.z[i] as f64));
            tree_args.push(&JsValue::from(self.species[i]));
            tree_args.push(&JsValue::from_f64(self.height[i] as f64));
            tree_args.push(&JsValue::from(self.status[i]));
            tree_args.push(&JsValue::from(self.stand_number[i]));

            let this = JsValue::null();
            callback.call1(&this, &tree_args).unwrap();
        }
    }

    // Getters for each field returning JavaScript arrays
    #[wasm_bindgen]
    pub fn x(&self) -> Array {
        let array = Array::new();
        for &val in &self.x {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn y(&self) -> Array {
        let array = Array::new();
        for &val in &self.y {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn z(&self) -> Array {
        let array = Array::new();
        for &val in &self.z {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn species(&self) -> Array {
        let array = Array::new();
        for &val in &self.species {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn height(&self) -> Array {
        let array = Array::new();
        for &val in &self.height {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn status(&self) -> Array {
        let array = Array::new();
        for &val in &self.status{
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }

    #[wasm_bindgen]
    pub fn stand_number(&self) -> Array {
        let array = Array::new();
        for &val in &self.stand_number {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }
}
