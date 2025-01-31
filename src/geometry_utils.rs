use crate::forest_property::compartment::find_stands_in_bounding_box;
use crate::forest_property::forest_property_data::{ForestPropertyData, ForestPropertyDataSchema};
use crate::forest_property::stand::Stand;
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::forest_property::tree::Tree;
use crate::jittered_hexagonal_sampling::{GridOptions, JitteredHexagonalGridSampling};
use crate::projection::{Projection, CRS};
use crate::shared_buffer::SharedBuffer;

use geo_types::Polygon;
use geo::{coord, Area, BooleanOps, BoundingRect, Coord, LineString};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use wasm_bindgen::JsValue;
use wasm_bindgen_futures::wasm_bindgen;
use wasm_bindgen::prelude::wasm_bindgen;
use web_sys::console::log_1;
use core::f32::consts::PI;

// Get minimum and maximum x and y coordinates of a polygon
pub fn get_min_max_coordinates(p: &Polygon<f64>) -> (f64, f64, f64, f64) {
    let rect = p.bounding_rect().unwrap();
    
    let min_x = rect.min().x;
    let max_x = rect.max().x;
    let min_y = rect.min().y;
    let max_y = rect.max().y;

    (min_x, max_x, min_y, max_y)
}

// Get the bounding box of the whole map
pub fn get_coords_of_map(xml: &str) -> (f64, f64, f64, f64) {
    let property = ForestPropertyData::from_xml_str(xml);
    let mut all_stands = property.real_estates.unwrap().real_estate[0].get_stands();

    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;

    for stand in all_stands.iter_mut() {
        let polygon = stand.computed_polygon.to_owned().unwrap();
        let (p_min_x, p_max_x, p_min_y, p_max_y) = get_min_max_coordinates(&polygon);

        if p_min_x < min_x {
            min_x = p_min_x;
        }
        if p_max_x > max_x {
            max_x = p_max_x;
        }
        if p_min_y < min_y {
            min_y = p_min_y;
        }
        if p_max_y > max_y {
            max_y = p_max_y;
        }
    }
    
    (min_x, max_x, min_y, max_y)
}

pub fn polygon_to_wgs84(p: &Polygon) -> Polygon {
    let proj = Projection::new(CRS::Epsg3067, CRS::Epsg4326);
    let mut coords: Vec<Coord<f64>> = Vec::new();

    for coordinate in p.exterior().points() {
        let e: f64 = coordinate.x();
        let n: f64 = coordinate.y();
        
        let (lon, lat) = proj.transform(e, n);
        coords.push(Coord { x: lon, y: lat });
    }

    Polygon::new(LineString::from(coords), vec![])
}

pub fn generate_radius(
    total_stem_count: u32, 
    area: f32
) -> f32 {

    let total_trees = total_stem_count as f32 * area / 10000.0;

    let mut ratio_fix = 1.3;

    if total_trees < 250.0 {
        ratio_fix = ((total_trees / 250.0) * 0.6) + 1.3;
    }
    let square_to_circle_ratio = 1.273 / ratio_fix;

    let tree_needed_area = area / total_trees / square_to_circle_ratio;

    // Calculate the radius based on the mean height of the tree species
    (tree_needed_area / PI).sqrt()
}


// Get the area ratio of the stand's clipped polygon in the bounding box to the original polygon
#[wasm_bindgen]
pub fn get_area_ratio(
    xml_content: &str,
    stand_id: u16,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
) -> JsValue {
    let bbox = Polygon::new(
        LineString(vec![
            coord!(x: x_min, y: y_min),
            coord!(x: x_max, y: y_min),
            coord!(x: x_max, y: y_max),
            coord!(x: x_min, y: y_max),
            coord!(x: x_min, y: y_min),
        ]),
        vec![],
    );

    let property = ForestPropertyData::from_xml_str(xml_content);
    let real_estate = property.real_estates.unwrap().real_estate[0].clone();
    let all_stands = real_estate.get_stands();

    let stands = find_stands_in_bounding_box(&all_stands, &bbox).unwrap();

    let stand = stands
        .iter()
        .find(|stand| stand.id.parse::<f64>().unwrap() == stand_id as f64)
        .unwrap();

    let polygon = stand.computed_polygon.to_owned().unwrap();
    let original_area = polygon.unsigned_area();
    let clipped_polygon = polygon.intersection(&bbox);
    let clipped_area = clipped_polygon.unsigned_area();

    let area_ratio = clipped_area / original_area;
    log_1(&format!("Area ratio: {}", area_ratio).into());
    JsValue::from(area_ratio)
}

// Generates random trees in stand for all strata with jittered grid sampling
pub fn generate_random_trees(
    p: &Polygon, 
    strata: &TreeStrata, 
    area_ratio: f64, 
    stand_id: f64
) -> Vec<Tree> {

    let total_stem_count = strata.tree_stratum.iter().fold(0, |mut acc: u32, f| {
        acc += f.stem_count;
        acc
    });

    let trees = strata
        .tree_stratum
        .par_iter()
        .map(|stratum| {
            let tree_amount = (stratum.stem_count as f64) * area_ratio;
            let amount = tree_amount.round() as u32;

            let mut radius = generate_radius(
                total_stem_count,
                stratum.basal_area
            );
            
            radius *= 0.00001;

            // Jittered Grid Version 2
            let rng = rand::thread_rng();
            let options = GridOptions {
                polygon: p.to_owned(),
                radius: (radius).into(),
                jitter: Some(0.6666),
                point_limit: Some(amount as usize),
            };

            let mut grid = JitteredHexagonalGridSampling::new(rng, options);

            let points =  grid.fill();
            
            if points.len() != 0 && points.len() < amount as usize {
                println!("Generated {} / {} trees for stratum with basal area {}, stem count {}, mean height {}.", points.len(), amount, stratum.basal_area, stratum.stem_count, stratum.mean_height);
            }

            let trees_strata: Vec<Tree> = points.iter().map(|pair: &[f64; 2]| {
                Tree::new(
                    stand_id,
                    stratum.tree_species,
                    stratum.mean_height,
                    (pair[0], pair[1], 0.0),
                    None
                )
            }).collect();
            
            trees_strata
        })
        .flatten();

    trees.collect()
}

pub fn generate_random_trees_no_strata(
    p: &Polygon,
    stand: &Stand,
    area_ratio: f64, 
    stand_id: f64
) -> Vec<Tree> {

    let stem_count = stand.summary_stem_count().unwrap_or(0);
    let tree_species = stand.main_tree_species().unwrap_or(0);
    let basal_area = stand.summary_basal_area().unwrap_or(0.0);
    let mean_height = stand.summary_mean_height().unwrap_or(0.0);

    let tree_amount = (stem_count as f64) * area_ratio;
    let amount = tree_amount.round() as u32;

    let mut radius = generate_radius(
        stem_count,
        basal_area
    );
            
    radius *= 0.00001;

    // Jittered Grid Version 2
    let rng = rand::thread_rng();
    let options = GridOptions {
        polygon: p.to_owned(),
        radius: (radius).into(),
        jitter: Some(0.6666),
        point_limit: Some(amount as usize),
    };

    let mut grid = JitteredHexagonalGridSampling::new(rng, options);

    let points =  grid.fill();

    let trees: Vec<Tree> = points.iter().map(|pair: &[f64; 2]| {
        Tree::new(
            stand_id,
            tree_species,
            mean_height,
            (pair[0], pair[1], 0.0),
            None
        )
    }).collect();
    
    trees
}

// Generates random trees for all strata in stand with jittered grid sampling
pub fn generate_random_trees_into_buffer(
    stand_p: &Polygon,
    clipped_p: &Polygon,
    strata: &TreeStrata,
    stand_id: f64,
    buffer: &SharedBuffer, // Pass in the SharedBuffer to fill
    start_index: usize,
) -> usize {

    let total_stem_count = strata.tree_stratum.iter().fold(0, |mut acc: u32, f| {
        acc += f.stem_count;
        acc
    });

    let trees = strata
        .tree_stratum
        .par_iter()
        .map(|stratum| {

            // Calculate the area ratio of the clipped polygon to the original polygon
            let original_area = stand_p.unsigned_area();
            let clipped_area = clipped_p.unsigned_area();
            let area_ratio = clipped_area / original_area;

            let tree_amount = (stratum.stem_count as f64) * area_ratio;
            let amount = tree_amount.round() as u32;

            let mut radius = generate_radius(total_stem_count, stratum.basal_area);
            radius *= 0.00001;

            // Jittered Grid Version 2
            let rng = rand::thread_rng();
            let options = GridOptions {
                polygon: clipped_p.to_owned(),
                radius: (radius).into(),
                jitter: Some(0.6666),
                point_limit: Some(amount as usize),
            };

            let mut grid = JitteredHexagonalGridSampling::new(rng, options);
            let points = grid.fill();

            let trees_strata: Vec<Tree> = points
                .iter()
                .map(|pair: &[f64; 2]| {
                    Tree::new(
                        stand_id,
                        stratum.tree_species,
                        stratum.mean_height,
                        (pair[0], pair[1], 0.0),
                        None
                    )
                })
                .collect();

            trees_strata
        })
        .flatten()
        .collect::<Vec<Tree>>();

    // Insert the trees into the buffer
    for (i, tree) in trees.iter().enumerate() {
        let buffer_index = start_index + i;
        if i < buffer.len() / 6 {

            // Fill the buffer with tree data
            buffer.fill_tree(
                buffer_index,
                tree.stand_id(),
                tree.position().0,
                tree.position().1,
                tree.species(),
                tree.tree_height(),
            );
        } else {
            break; // Avoid overflowing the buffer
        }
    }

    trees.len() // Return the number of trees added to the buffer
}