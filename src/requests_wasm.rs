use crate::geometry_utils::{get_min_max_coordinates, generate_radius};
use crate::forest_property::forest_property_data::ForestPropertyData;
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::forest_property::tree::Tree;
use crate::forest_property::stand::Stand;
use crate::forest_property::compartment::{find_stands_in_bounding_box, CompartmentArea};
use crate::geojson_utils::all_compartment_areas_to_geojson;
use crate::shared_buffer::SharedBuffer;
use crate::jittered_hexagonal_sampling::{GridOptions, JitteredHexagonalGridSampling};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use geo::{coord, Area, LineString, Polygon, BooleanOps};
use geojson::{GeoJson, Value};
use reqwest_wasm::Client;
use reqwest::Error as ReqwestError;
use geojson::Error as GeoJsonError;
use std::fmt;
use std::error::Error;
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::JsValue;
use wasm_bindgen_futures::JsFuture;
use serde_json::Value as SerdeJsonValue;
use serde_wasm_bindgen;
use web_sys::console::log_1;
use serde::Serialize;
use std::slice::from_raw_parts;

#[derive(Debug)]
pub enum FetchError {
    Reqwest(ReqwestError),
    GeoJson(GeoJsonError),
}

impl fmt::Display for FetchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FetchError::Reqwest(err) => write!(f, "Reqwest error: {}", err),
            FetchError::GeoJson(err) => write!(f, "GeoJson error: {}", err),
        }
    }
}

impl std::error::Error for FetchError {}

impl From<ReqwestError> for FetchError {
    fn from(err: ReqwestError) -> Self {
        FetchError::Reqwest(err)
    }
}

impl From<GeoJsonError> for FetchError {
    fn from(err: GeoJsonError) -> Self {
        FetchError::GeoJson(err)
    }
}

pub fn geojson_to_polygons(geojson: &GeoJson) -> Vec<Polygon<f64>> {
    // Initialize a vector to store polygons
    let mut polygons = Vec::new();

    // Match on GeoJson to handle FeatureCollection
    if let GeoJson::FeatureCollection(collection) = geojson {
        for feature in &collection.features {
            // Ensure we are working with a valid Feature
            if let Some(geometry) = &feature.geometry {
                match &geometry.value {
                    Value::Polygon(polygon) => {
                        // Convert GeoJSON Polygon to geo crate Polygon
                        let exterior = polygon[0]
                            .iter()
                            .map(|point| (point[0], point[1]))
                            .collect::<Vec<_>>();
                        
                        // Create a geo crate Polygon
                        let poly = Polygon::new(LineString::from(exterior), vec![]);
                        polygons.push(poly);
                    }
                    _ => {
                        // Handle other geometry types if necessary
                        eprintln!("Skipping non-polygon geometry");
                    }
                }
            }
        }
    } 

    polygons
}

#[derive(Serialize)]
struct GeoJsonWithTreeCount {
    geojson: serde_json::Value,
    max_tree_count: usize,
    tree_count: usize,
    buffer_pointer: u64,
}

#[wasm_bindgen]
pub async fn geo_json_from_coords(
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
    xml_content: String,
) -> Result<JsValue, JsValue> {
    // Get the ForestPropertyData from the XML content
    let property = ForestPropertyData::from_xml_str(&xml_content);
    log_1(&"Got property".into());

    let mut bbox = Polygon::new(
        LineString(vec![
            coord!(x: min_x, y: min_y),
            coord!(x: max_x, y: min_y),
            coord!(x: max_x, y: max_y),
            coord!(x: min_x, y: max_y),
            coord!(x: min_x, y: min_y),
        ]),
        vec![],
    );

    let west = min_x;
    let south = min_y;
    let east = max_x;
    let north = max_y;

    let url_buildings = format!(
        "https://metne-test.onrender.com/geoserver/mml/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=mml:rakennus&maxFeatures=2000&outputFormat=application%2Fjson&BBOX={},{},{},{},EPSG:4326&srsName=EPSG:4326",
        west, south, east, north
    );
    let url_roads = format!(
        "https://metne-test.onrender.com/geoserver/mml/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=mml:tieviiva&bbox={},{},{},{},EPSG:4326&srsName=EPSG:4326&outputFormat=application/json",
        west, south, east, north
    );

    // Create HTTP client for async fetch
    let client = Client::new();

    // Fetch buildings GeoJSON
    let buildings_response = client.get(&url_buildings).send().await
        .map_err(|e| JsValue::from_str(&format!("Failed to fetch buildings: {}", e)))?;
    let buildings_text = buildings_response.text().await
        .map_err(|e| JsValue::from_str(&format!("Failed to read buildings response text: {}", e)))?;
    let buildings_geojson: GeoJson = serde_json::from_str(&buildings_text)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse buildings GeoJson: {}", e)))?;

    let buildings = geojson_to_polygons(&buildings_geojson);
    let buildings_count = buildings.len();
    log_1(&format!("Fetched {} buildings", buildings_count).into());

    // Exclude buildings from the bounding box
    for building in buildings.iter() {
        bbox = bbox.difference(building).0.first().unwrap().to_owned();
    }

    // Fetch roads GeoJSON
    let roads_response = client.get(&url_roads).send().await
        .map_err(|e| JsValue::from_str(&format!("Failed to fetch roads: {}", e)))?;
    let roads_text = roads_response.text().await
        .map_err(|e| JsValue::from_str(&format!("Failed to read roads response text: {}", e)))?;
    let roads_geojson: GeoJson = serde_json::from_str(&roads_text)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse roads GeoJson: {}", e)))?;

    // Get the ForestPropertyData and stands
    let real_estate = property.real_estates.real_estate[0].clone();
    let stands = real_estate.get_stands();

    // Get compartment areas in the bounding box and convert them to GeoJSON
    let compartment_areas = get_compartment_areas_in_bounding_box(stands, &bbox);
    let max_tree_count = compartment_areas.1;
    let tree_count = compartment_areas.2;
    let buffer_pointer = compartment_areas.3;
    let geojson = all_compartment_areas_to_geojson(compartment_areas.0, &buildings_geojson, &roads_geojson);
    log_1(&"Got geojson".into());

    // Create a combined struct with both the GeoJSON and tree_count
    let result = GeoJsonWithTreeCount {
        geojson: geojson.into(),
        max_tree_count,
        tree_count,
        buffer_pointer,
    };

    // Serialize the result to a JsValue to return to JavaScript
    let result_js_value = serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsValue::from_str(&format!("Failed to serialize result: {}", e)))?;

    Ok(result_js_value)
}

// Generates random trees for all strata with jittered grid sampling
pub fn generate_random_trees_into_buffer(
    p: &Polygon,
    strata: &TreeStrata,
    area_ratio: f64,
    buffer: &SharedBuffer, // Pass in the SharedBuffer to fill
    start_index: usize
) -> usize {
    let total_stem_count = strata.tree_stratum.iter().fold(0, |mut acc: u32, f| {
        acc += f.stem_count;
        acc
    });

    let mut tree_count = 0;

    let trees = strata
        .tree_stratum
        .par_iter()
        .map(|stratum| {
            let tree_amount = (stratum.stem_count as f64) * area_ratio;
            let amount = tree_amount.round() as u32;

            let mut radius = generate_radius(total_stem_count, stratum.basal_area);
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
            let points = grid.fill();

            if points.is_empty() {
                //println!("\tNo trees generated for stratum with basal area {}, stem count {}, mean height {}", stratum.basal_area, stratum.stem_count, stratum.mean_height);
            } else if points.len() < amount as usize {
                println!(
                    "Generated {} / {} trees for stratum with basal area {}, stem count {}, mean height {}.",
                    points.len(), amount, stratum.basal_area, stratum.stem_count, stratum.mean_height
                );
            }

            let trees_strata: Vec<Tree> = points
                .iter()
                .map(|pair: &[f64; 2]| {
                    Tree::new(stratum.tree_species, stratum.mean_height, (pair[0], pair[1], 0.0))
                })
                .collect();

            trees_strata
        })
        .flatten()
        .collect::<Vec<Tree>>();
 
    // Insert the trees into the buffer
    for (i, tree) in trees.iter().enumerate() {
        let buffer_index = start_index + i;
        if i < buffer.len() / 3 {
            // Fill the buffer with x, y, and species
            buffer.fill_tree(buffer_index, tree.position().0, tree.position().1, tree.species());
            tree_count += 1;
        } else {
            break; // Avoid overflowing the buffer
        }
    }

    tree_count // Return the number of trees added to the buffer
}

// Get compartment areas in a bounding box.
pub fn get_compartment_areas_in_bounding_box(
    all_stands: Vec<Stand>,
    bbox: &Polygon,
) -> (Vec<CompartmentArea>, usize, usize, u64) {
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

    // Create a shared buffer to store the generated trees
    let buffer = SharedBuffer::new(max_tree_count as usize);

    // If there are stands in the bounding box, generate random trees for each stand
    if let Some(stands) = stands {
        let mut compartment_areas = Vec::new();
        let mut total_tree_count = 0;

        let mut buffer_index = 0;
        for stand in stands {
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

            // Generate trees and save them to the buffer if strata exist
            let mut tree_count = 0;
            if let Some(strata) = strata {
                tree_count = generate_random_trees_into_buffer(&clipped_polygon, &strata, area_ratio, &buffer, buffer_index);
                buffer_index += tree_count;
                log_1(&format!("Generated {} trees for stand {}", tree_count, stand.stand_basic_data.stand_number).into());
            }
            total_tree_count += tree_count;

            // Add to the compartment areas list
            compartment_areas.push(CompartmentArea {
                stand_number: stand.stand_basic_data.stand_number.to_string(),
                polygon: clipped_polygon,
            });
        }

        // Get a slice of the buffer
        let buffer_slice: &[f64] = unsafe {
            std::slice::from_raw_parts(buffer.ptr(), buffer.len())
        };  

        log_1(&"Buffer contains:".into());
        for (i, value) in buffer_slice.iter().enumerate() {
            if i % 3 == 0 && buffer_slice[i + 2] != 0.0 {
                let buffer_info = format!("Tree {}: x = {}, y = {}, species = {}", i / 3, buffer_slice[i], buffer_slice[i + 1], buffer_slice[i + 2]);
                log_1(&buffer_info.into());
            }
        }
        
        let mut buffer_pointer_string = format!("{:p}", buffer.ptr());
        buffer_pointer_string = (&buffer_pointer_string[2..]).to_string();

        log_1(&format!("Hexadecimal Buffer pointer in rust: {}", buffer_pointer_string).into());
        let buffer_pointer = hexadecimal_to_decimal(&buffer_pointer_string).unwrap();
        log_1(&format!("Decimal Buffer pointer in rust: {}", buffer_pointer).into());
        (compartment_areas, max_tree_count as usize, total_tree_count, buffer_pointer)
    } else {
        (vec![], 0, 0, 0)
    }
}

pub fn hexadecimal_to_decimal(hexadecimal_str: &str) -> Result<u64, &'static str> {
    if hexadecimal_str.is_empty() {
        return Err("Empty input");
    }

    for hexadecimal_str in hexadecimal_str.chars() {
        if !hexadecimal_str.is_ascii_hexdigit() {
            return Err("Input was not a hexadecimal number");
        }
    }

    match u64::from_str_radix(hexadecimal_str, 16) {
        Ok(decimal) => Ok(decimal),
        Err(_e) => Err("Failed to convert to hexadecimal"),
    }
}