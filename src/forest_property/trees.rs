use crate::{forest_property::tree::Tree, geojson_utils::convert_tree_to_feature};

use geo::Point;
use geojson::{Feature, FeatureCollection, GeoJson};
use serde_wasm_bindgen::to_value;
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
    stand_id: Vec<u16>,
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
        stand_id: Vec<u16>,
    ) -> Trees {
        Trees {
            x,
            y,
            z,
            species,
            height,
            status,
            stand_id,
        }
    }

    pub fn insert(&mut self, tree: Tree) {
        self.x.push(tree.position().0 as f32);
        self.y.push(tree.position().1 as f32);
        self.z.push(tree.position().2 as f32);
        self.species.push(tree.species());
        self.height.push(tree.tree_height());
        self.status.push(tree.tree_status() as u8);
        self.stand_id.push(tree.stand_id() as u16);
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
            tree_args.push(&JsValue::from(self.stand_id[i]));

            let this = JsValue::null();
            callback.call1(&this, &tree_args).unwrap();
        }
    }

    pub fn to_geojson(&self) -> JsValue {
        let mut trees = Vec::new();
        let trees_length = self.x.len();

        for i in 0..trees_length {
            let stand_id = self.stand_id[i] as f64;
            let species = self.species[i];
            let height = self.height[i];
            let position = (self.x[i] as f64, self.y[i] as f64, self.z[i] as f64);
            let status = self.status[i] as f64;

            let tree = Tree::new(stand_id, species, height, position, Some(status));
            trees.push(tree);
        }

        // Convert the trees to GeoJSON
        let tree_features: Vec<Feature> = trees.iter().map(|tree| convert_tree_to_feature(tree)).collect();

        // Create the GeoJSON FeatureCollection
        let feature_collection = FeatureCollection {
            features: tree_features,
            bbox: None,
            foreign_members: None,
        };

        // Return a GeoJson object
        let geojson = GeoJson::FeatureCollection(feature_collection);

        let geojson_string = geojson.to_string();
        
        // Convert the String to JsValue
        JsValue::from_str(&geojson_string)
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
    pub fn stand_id(&self) -> Array {
        let array = Array::new();
        for &val in &self.stand_id {
            array.push(&JsValue::from_f64(val as f64));
        }
        array
    }
}
