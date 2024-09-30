use crate::forest_property::tree::Tree;
use crate::geometry_utils::generate_random_trees;
use super::stand::Stand;

use geo::{Polygon, Area, BooleanOps};
use geo::Intersects;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

// Struct that represents a stand of trees
#[derive(Debug, Clone)]
pub struct Compartment {
    pub stand_number: String,
    pub trees: Vec<Tree>,
    pub polygon: Polygon,
}

impl Compartment {
    pub fn new(stand_number: String, trees: Vec<Tree>, polygon: Polygon) -> Self {
        Compartment {
            stand_number,
            trees,
            polygon,
        }
    }

    pub fn stand_number(&self) -> &String {
        &self.stand_number
    }

    pub fn trees(&self) -> &Vec<Tree> {
        &self.trees
    }

    pub fn polygon(&self) -> &Polygon {
        &self.polygon
    }

    // Polygon clipping to bounding box
    pub fn clip_polygon_to_bounding_box(&self, bbox: &Polygon) -> Option<Polygon> {
        let clipped = self.polygon.intersection(bbox);

        if clipped.0.is_empty() {
            println!("Polygon is empty");
            None
        } else {
            let p = clipped.0.first().unwrap();
            Some(p.to_owned())
        }
    }

    // Get trees in a bounding box
    pub fn trees_in_bounding_box(&self, min_x: f64, max_x: f64, min_y: f64, max_y: f64) -> Vec<&Tree> {
        self.trees.iter().filter(|tree| {
            let (x, y, _) = tree.position();
            x >= min_x && x <= max_x && y >= min_y && y <= max_y  // Keep the tree if it is inside the bounding box
        }).collect()
    }
}

pub fn find_stands_in_bounding_box<'a>(stands: &'a Vec<Stand>, bbox: &'a Polygon) -> Option<Vec<&'a Stand>> {

    // Collect the stands that intersect with the bounding box
    let intersecting_stands: Vec<&Stand> = stands.iter().filter(|stand| {
        let (exterior, _) = stand.get_geometries();
        bbox.intersects(&exterior)
    }).collect();  // Collect the stands that intersect with the bounding box

    if intersecting_stands.is_empty() {
        println!("No stands found in the bounding box");
        None
    } else {
        Some(intersecting_stands)
    }
}

// Get compartments in a bounding box.
pub fn get_compartments_in_bounding_box(
    all_stands: Vec<Stand>,
    bbox: &Polygon
) -> Vec<Compartment> {
    // Find stands in the bounding box
    let stands = find_stands_in_bounding_box(&all_stands, bbox);

    // If there are stands in the bounding box, generate random trees for each stand
    if !&stands.is_none() {
        let compartments: Vec<Compartment> = stands.unwrap()
            .into_par_iter()
            .map(|stand| {
                let polygon = stand.computed_polygon.to_owned().unwrap();
                let strata = stand.get_strata();

                // Clip the stand's polygon to the bounding box
                let intersected_polygons = polygon.intersection(bbox).0;
                let clipped_polygon = intersected_polygons.first()
                    .expect("Intersection result should contain at least one polygon")
                    .to_owned();
                
                // Calculate the area ratio of the clipped polygon to the original polygon
                let original_area = polygon.unsigned_area();
                let clipped_area = clipped_polygon.unsigned_area();
                let area_ratio = clipped_area / original_area;

                // Generate trees if strata exist
                let trees = if let Some(strata) = strata {
                    generate_random_trees(&clipped_polygon, &strata, area_ratio)
                } else {
                    vec![]
                };

                // Create and return the compartment
                Compartment {
                    stand_number: stand.stand_basic_data.stand_number.to_string(),
                    trees,
                    polygon: clipped_polygon,
                }
            })
            .collect();

        compartments
    } else {
        vec![]
    }
}

pub struct CompartmentArea {
    pub stand_number: String,
    pub polygon: Polygon,
}

