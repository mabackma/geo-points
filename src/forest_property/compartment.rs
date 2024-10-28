use crate::forest_property::tree::Tree;
use crate::geometry_utils::generate_random_trees;
use crate::projection::{Projection, CRS};
use crate::requests_wasm::{check_simulation, get_cutting_volume, get_simulation_strata, OperationType};
use super::stand::Stand;

use geo::{Area, BooleanOps, Contains, Coord, Polygon};
use geo::Intersects;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use web_sys::console::log_1;

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

    pub fn operate_compartment(&mut self, op_type: OperationType, area_polygons: Vec<Polygon>) {
        // Make a projection from epsg3067 to epsg4326
        let projection = Projection::new(CRS::Epsg3067, CRS::Epsg4326);
        let areas = projection.polygons_3067_to_4326(area_polygons);

        let cutting_volume = get_cutting_volume(&op_type);

        if op_type == OperationType::Cutting(cutting_volume) {
            log_1(&format!("Cutting operation with volume: {}", cutting_volume).into());
            log_1(&format!("Areas: {:?}", areas).into());
            // Count the trees within the specified areas
            let tree_count = self.trees.iter().filter(|tree| {
                let (x, y, _) = tree.position();
                areas.iter().any(|area| area.contains(&geo::Point::new(x, y)))
            }).count();

            // Calculate the number of trees to cut
            let trees_to_cut = ((cutting_volume / 100.0) * tree_count as f64).round() as usize;

            log_1(&format!("Cutting {} / {} trees", trees_to_cut, tree_count).into());
            let mut cut_count = 0;
            for tree in self.trees.iter_mut() {

                // Check if the tree is within any of the specified areas
                let (x, y, _) = tree.position();
                if areas.iter().any(|area| area.contains(&geo::Point::new(x, y))) {

                    // Cut the tree if the number of trees to cut has not been reached
                    if cut_count < trees_to_cut {
                        tree.cut_tree();
                        cut_count += 1;
                    } else {
                        break; 
                    }
                }
            }
        } 
        
        if op_type == OperationType::Thinning(cutting_volume) {
            // Count the trees within the specified areas
            let tree_count = self.trees.iter().filter(|tree| {
                let (x, y, _) = tree.position();
                areas.iter().any(|area| area.contains(&geo::Point::new(x, y)))
            }).count();

            // implement thinning operation
        } 

        if check_simulation(&op_type) {
            // implement simulation operation
            let strata = get_simulation_strata(&op_type);
        }
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
                let stand_number = stand.stand_basic_data.stand_number as f64;

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
                    generate_random_trees(&clipped_polygon, &strata, area_ratio, stand_number)
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

// Get compartment areas in a bounding box.
// Returns a tuple of compartment areas, max tree count, tree count, and buffer pointer in decimal
pub fn get_compartment_areas_in_bounding_box(
    all_stands: Vec<Stand>,
    bbox: &Polygon,
) -> Option<(Vec<CompartmentArea>, usize, usize)> {
    // Find stands in the bounding box
    let stands = find_stands_in_bounding_box(&all_stands, bbox);

    // Count the total number of trees in the bounding box
    let mut max_tree_count = 0;
    if let Some(stands) = &stands {
        for stand in stands {
            let strata = stand.get_strata();

            if let Some(strata) = strata {
                let strata_stem_count = strata.tree_stratum.iter().fold(0, |mut acc: u32, f| {
                    acc += f.stem_count;
                    acc
                });
                max_tree_count += strata_stem_count;
            }
        }
    }

    // If there are stands in the bounding box, generate random trees for each stand
    if let Some(stands) = stands {
        let mut compartment_areas = Vec::new();
        let mut total_tree_count = 0;

        for stand in stands {
            let polygon = stand.computed_polygon.to_owned().unwrap();

            // Clip the stand's polygon to the bounding box
            let intersected_polygons = polygon.intersection(bbox).0;
            let clipped_polygon = intersected_polygons
                .first()
                .expect("Intersection result should contain at least one polygon")
                .to_owned();

            total_tree_count += stand.summary_stem_count().unwrap_or(0) as usize;
            // Add to the compartment areas list
            compartment_areas.push(CompartmentArea {
                stand_number: stand.stand_basic_data.stand_number.to_string(),
                polygon: clipped_polygon,
            });
        }
        return Some((compartment_areas, max_tree_count as usize, total_tree_count));
    } else {
        None
    }
}