use std::collections::HashSet;

use crate::forest_property::tree::Tree;
use crate::geometry_utils::generate_random_trees;
use crate::projection::{Projection, CRS};
use crate::requests_wasm::OperationType;
use super::stand::Stand;
use super::tree_stand_data::TreeStrata;

use geo::{BooleanOps, Contains, HaversineDistance, Point, Polygon};
use geo::algorithm::area::Area;
use geo::Intersects;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use web_sys::console::log_1;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone)]
struct HashablePoint(Point<f64>);

impl Hash for HashablePoint {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Access x and y coordinates directly from Point
        let (x, y) = (self.0.x(), self.0.y());
        x.to_bits().hash(state); // Hash the x-coordinate
        y.to_bits().hash(state); // Hash the y-coordinate
    }
}

impl PartialEq for HashablePoint {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Eq for HashablePoint {}

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

        let cutting_volume = op_type.get_cutting_volume();

        // Count the trees within the specified areas
        let tree_count = self.trees.iter().filter(|tree| {
            let (x, y, _) = tree.position();
            areas.iter().any(|area| area.contains(&geo::Point::new(x, y)))
        }).count();
    
        // Calculate the number of trees to cut
        let trees_to_cut = ((cutting_volume / 100.0) * tree_count as f64).round() as usize;
        
        log_1(&format!("{} operation with volume {} on {} / {} trees", op_type.get_type(), cutting_volume, trees_to_cut, tree_count).into());
        log_1(&format!("Areas: {:?}", areas).into());

        if op_type == OperationType::Cutting(cutting_volume) {
            self.cutting_operation(trees_to_cut, &areas);
        } 
        
        if op_type == OperationType::Thinning(cutting_volume) {
            let trees_left = tree_count - trees_to_cut;
            let thinning_distance = self.calculate_thinning_distance(&areas, trees_left as f64);

            self.thinning_operation(thinning_distance, trees_to_cut, &areas);
        } 

        if op_type.check_simulation() {
            let strata = op_type.get_simulation_strata();
            self.simulation(strata);
        }
    }

    fn calculate_thinning_distance(&self, areas: &Vec<Polygon>, trees_left: f64) -> f64 {
        let mut radius = 0.1; 

        // Project areas to EPSG:3067 for area calculation in m^2
        let proj = Projection::new(CRS::Epsg4326, CRS::Epsg3067);
        let areas = proj.polygons_4326_to_3067(areas.clone());
        let total_area = areas.iter().fold(0.0, |acc, area| acc + area.unsigned_area());
        
        let ratio_fix = 2.0;
        let square_to_circle_ratio = 1.273 / ratio_fix;

        let tree_needed_area = total_area / trees_left / square_to_circle_ratio;
        log_1(&format!("Trees left after thinning: {}", trees_left).into());
        log_1(&format!("Total area: {}", total_area).into());
        log_1(&format!("Tree-needed area: {}", tree_needed_area).into());
        radius = (tree_needed_area / std::f64::consts::PI).sqrt();
        radius = radius.max(0.1); 

        radius * 2.0
    }
    
    fn cutting_operation(&mut self, trees_to_cut: usize, areas: &Vec<Polygon>) {
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

    // Thinning Logic: The trees that are added to all_removed_positions will be the ones to be removed
    // Points to be removed have neighbors that are too close together.
    fn thinning_operation(&mut self, thinning_distance: f64, trees_to_cut: usize, areas: &Vec<Polygon<f64>>) {
        log_1(&format!("Thinning operation with distance: {} and {} trees to cut", thinning_distance, trees_to_cut).into());
        let mut all_removed_positions = HashSet::new();

        for area in areas {
            let trees_in_area: Vec<_> = self.trees.iter()
                .filter(|tree| {
                    let (x, y, _) = tree.position();
                    area.contains(&Point::new(x, y))
                })
                .cloned()
                .collect();

            let mut removed_positions_in_area = HashSet::new();

            for (i, tree) in trees_in_area.iter().enumerate() {
                let (tree_x, tree_y, _) = tree.position();
                let tree_point = HashablePoint(Point::new(tree_x, tree_y)); // Use the wrapper type

                // Skip trees that are already marked for removal in this area
                if removed_positions_in_area.contains(&tree_point) {
                    continue;
                }

                // Check all remaining trees in the area for thinning distance
                for other_tree in trees_in_area.iter().skip(i + 1) {
                    let (other_x, other_y, _) = other_tree.position();
                    let other_point = HashablePoint(Point::new(other_x, other_y)); // Use the wrapper type

                    // Stop thinning if we exceed the trees_to_cut limit
                    if all_removed_positions.len() + removed_positions_in_area.len() >= trees_to_cut {
                        break;
                    }

                    if HaversineDistance::haversine_distance(&tree_point.0, &other_point.0) < thinning_distance {
                        removed_positions_in_area.insert(tree_point);
                        break;
                    }
                }
            }

            all_removed_positions.extend(removed_positions_in_area);

            // Stop thinning if we exceed the trees_to_cut limit
            if all_removed_positions.len() >= trees_to_cut {
                break; 
            }
        }

        log_1(&format!("Removing {} trees", all_removed_positions.len()).into());

        // Cut the trees based on `all_removed_positions`
        for tree in self.trees.iter_mut() {
            let (tree_x, tree_y, _) = tree.position();
            let tree_point = HashablePoint(Point::new(tree_x, tree_y));
            
            if all_removed_positions.contains(&tree_point) {
                tree.cut_tree();
            }
        }
    }

    fn simulation(&mut self, strata: TreeStrata) {
        log_1(&format!("Simulation operation with strata: {:#?}", strata).into());
        log_1(&format!("Trees in the compartment: {}", self.trees.len()).into());
        self.trees = generate_random_trees(&self.polygon, &strata, 1.0, self.stand_number.parse().unwrap());
        log_1(&format!("Trees in the compartment after simulation: {}", self.trees.len()).into());
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