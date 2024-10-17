use crate::forest_property::compartment::{
    find_stands_in_bounding_box, get_compartments_in_bounding_box, Compartment, CompartmentArea,
};
use crate::forest_property::forest_property_data::{ForestPropertyData, RealEstate, TreeStratum};
use crate::forest_property::stand::{self, Stand};
use crate::forest_property::tree::Tree;
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::geojson_utils::all_compartment_areas_to_geojson;
use crate::geometry_utils::{generate_radius, get_min_max_coordinates};
use crate::jittered_hexagonal_sampling::{GridOptions, JitteredHexagonalGridSampling};
use crate::shared_buffer::SharedBuffer;

use geo::{coord, point, Area, BooleanOps, Contains, Line, LineString, Polygon};
use geojson::{Error as GeoJsonError, JsonObject};
use geojson::{GeoJson, Value};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use reqwest::Error as ReqwestError;
use reqwest_wasm::Client;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value as SerdeJsonValue};
use serde_wasm_bindgen;
use std::error::Error;
use std::fmt;
use std::slice::from_raw_parts;
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
use web_sys::console::log_1;
use web_sys::js_sys::{Float32Array, JsString};

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

// Struct to hold GeoJSON data, max tree count, tree count, and buffer pointer
#[derive(Serialize)]
struct GeoJsonWithTreeCount {
    geojson: serde_json::Value,
    max_tree_count: usize,
    tree_count: usize,
    buffer_pointer: u64,
}
#[derive(Serialize,Deserialize, Clone, PartialEq)]
enum OperationType {
    Thinning,
    Cutting,
    Simulation,
}

#[derive(Serialize,Deserialize, Clone, PartialEq)]
#[wasm_bindgen]
struct Operation {
    operation_type: OperationType,
    cutting_areas: Vec<Polygon>,
}

#[derive(Serialize,Deserialize, Clone, PartialEq)]
#[wasm_bindgen]
struct StandOperation {
    id: u32,
    stand_id: u32,
    operations: Vec<TreeStrata>,
    active_operation: u32,
}

#[wasm_bindgen]
#[derive(Serialize,Deserialize, Clone, PartialEq)]
pub struct VirtualForest {
    property: ForestPropertyData,
    selected_realestate: u32,
    stand_operations: Vec<StandOperation>,
    retention_zones: Vec<Polygon>,
    roads: Option<Vec<LineString>>,
    buildings: Option<Vec<Polygon>>,
    water: Option<Vec<Polygon>>
}

#[wasm_bindgen]
#[derive(Serialize,Deserialize, Debug, Clone, PartialEq)]
pub struct RealEstateValue {
    pub id: u32,
    pub index: usize,
    pub group_number: u16,
    pub area_number: u16,
}

#[wasm_bindgen]
impl VirtualForest {
    #[wasm_bindgen(constructor)]
    pub fn new(xml: &str) -> Self {
        let property = ForestPropertyData::from_xml_str(xml);

        VirtualForest {
            property,
            selected_realestate: 0,
            stand_operations: vec![],
            retention_zones: vec![],
            roads: None,
            water: None,
            buildings: None
        }
    }

    #[wasm_bindgen]
    pub fn to_json(&self) -> JsValue {
   
        let data = serde_json::to_string(&self).unwrap();
        /*         log_1(&format!("property serialized {}", data).into()); */
        JsValue::from(data)
        // JsString::from(JsValue::self.property)
        // self.property.to
    }

    #[wasm_bindgen]
    pub fn from_json(json: &str) -> Self {


        let property : Self = serde_json::from_str(json).expect("Error parsing json file");


        property


    }

    #[wasm_bindgen]
    pub fn get_realestates(&self) -> Vec<RealEstateValue> {
        let realestates = self
            .property
            .real_estates
            .real_estate
            .iter()
            .enumerate()
            .map(|(i, f)| RealEstateValue {
                index: i,
                id: f.id,
                area_number: f.area_number,
                group_number: f.group_number,
            })
            .collect::<Vec<RealEstateValue>>();
        realestates
        /* JsValue::from(serde_json::to_string(&realestates).unwrap()) */
    }

    #[wasm_bindgen]
    pub fn set_realestate(&mut self, index: u32) {
        self.selected_realestate = index;
    }

    #[wasm_bindgen]
    pub fn get_selected_realestate(&self) -> Result<RealEstateValue, JsValue> {
        if let Some((index, result)) = self
            .property
            .real_estates
            .real_estate
            .iter()
            .enumerate()
            .find(|(i, _f)| *i == self.selected_realestate as usize)
        {
            Ok(RealEstateValue {
                id: result.id,
                index,
                group_number: result.group_number,
                area_number: result.area_number,
            })
        } else {
            Err(JsValue::NULL)
        }
    }

    fn _get_selected_realestate(&self) -> Option<RealEstate> {
        if let Some((index, result)) = self
            .property
            .real_estates
            .real_estate
            .iter()
            .enumerate()
            .find(|(i, _f)| *i == self.selected_realestate as usize)
        {
            Some(result.to_owned())
        } else {
            None
        }
    }

    #[wasm_bindgen]
    pub fn get_stand_by_id(&self, id: String) -> Result<JsValue, JsValue> {
        if let Some(stand) = self
            ._get_selected_realestate()
            .unwrap()
            .get_stands()
            .iter()
            .find(|f| f.id == id)
        {
            return Ok(JsValue::from(
                serde_json::to_string(&stand.to_owned()).expect("Stand should be serializable"),
            ));
        } else {
            Err(JsValue::NULL)
        }
    }

    #[wasm_bindgen]
    pub fn generate_trees_bbox(&self, min_x: f64, max_x: f64, min_y: f64, max_y: f64) -> Vec<f32> {
        let bbox = Polygon::new(
            LineString(vec![
                coord!(x: min_x, y: min_y),
                coord!(x: max_x, y: min_y),
                coord!(x: max_x, y: max_y),
                coord!(x: min_x, y: max_y),
                coord!(x: min_x, y: min_y),
            ]),
            vec![],
        );

        let stands = self._get_selected_realestate().unwrap().get_stands();

        let compartments = get_compartments_in_bounding_box(stands, &bbox);

        let trees = compartments
            .iter()
            .enumerate()
            .flat_map(|(i, f)| {
                f.trees
                    .iter()
                    .flat_map(|tree| {
                        let (x, y, z) = tree.position();
                        let specie = tree.species();
                        let height = tree.tree_height();
                        let status = tree.tree_status();
                        let stand_number = tree.stand_number();
                        vec![
                            x as f32,
                            y as f32,
                            z as f32,
                            specie as f32,
                            height,
                            status as f32,
                            stand_number as f32,
                        ]
                    })
                    .collect::<Vec<f32>>()
            })
            .collect::<Vec<f32>>();

        return trees;
    }

    // Fetches GeoJSON data from the given bounding box and XML content
    // Returns a GeoJsonWithTreeCount struct as a JsValue
    #[wasm_bindgen]
    pub async fn geo_json_from_coords(
        &self,
        min_x: f64,
        max_x: f64,
        min_y: f64,
        max_y: f64,
    ) -> Result<JsValue, JsValue> {
        /*         // Get the ForestPropertyData from the XML content
        let property = ForestPropertyData::from_xml_str(xml_content);

        log_1(&"Got property".into()); */

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
        let buildings_response = client
            .get(&url_buildings)
            .send()
            .await
            .map_err(|e| JsValue::from_str(&format!("Failed to fetch buildings: {}", e)))?;
        let buildings_text = buildings_response.text().await.map_err(|e| {
            JsValue::from_str(&format!("Failed to read buildings response text: {}", e))
        })?;
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
        let roads_response = client
            .get(&url_roads)
            .send()
            .await
            .map_err(|e| JsValue::from_str(&format!("Failed to fetch roads: {}", e)))?;
        let roads_text = roads_response.text().await.map_err(|e| {
            JsValue::from_str(&format!("Failed to read roads response text: {}", e))
        })?;
        let roads_geojson: GeoJson = serde_json::from_str(&roads_text)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse roads GeoJson: {}", e)))?;

        // Get the ForestPropertyData and stands
        let stands = self._get_selected_realestate().unwrap().get_stands();
        // let stands = real_estate.get_stands();

        // Get compartment areas in the bounding box and convert them to GeoJSON
        if let Some(compartment_areas) = get_compartment_areas_in_bounding_box(stands, &bbox) {
            /*         let max_tree_count = compartment_areas.1;
            let tree_count = compartment_areas.2;
            let buffer_pointer = compartment_areas.3; */
            let geojson = all_compartment_areas_to_geojson(
                compartment_areas.0,
                &buildings_geojson,
                &roads_geojson,
            );
            /*   log_1(&"Got geojson".into()); */

            // Create a combined struct with both the GeoJSON and tree_count
            /*        let result = GeoJsonWithTreeCount {
                       geojson: geojson.into(),
                       max_tree_count,
                       tree_count,
                       buffer_pointer,
                   };
            */
            // Serialize the result to a JsValue to return to JavaScript
            /*       let result_js_value = serde_wasm_bindgen::to_value(&result)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize result: {}", e)))?; */

            Ok(JsValue::from(geojson.to_string()))
        } else {
            Err(JsValue::from_str("could not create compartments"))
        }
    }
}

// Get the area ratio of the stand's clipped polygon in the bounding box to the original polygon
#[wasm_bindgen]
pub fn get_area_ratio(
    xml_content: &str,
    stand_number: u16,
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
    let real_estate = property.real_estates.real_estate[0].clone();
    let all_stands = real_estate.get_stands();

    let stands = find_stands_in_bounding_box(&all_stands, &bbox).unwrap();

    let stand = stands
        .iter()
        .find(|stand| stand.stand_basic_data.stand_number == stand_number)
        .unwrap();

    let polygon = stand.computed_polygon.to_owned().unwrap();
    let original_area = polygon.unsigned_area();
    let clipped_polygon = polygon.intersection(&bbox);
    let clipped_area = clipped_polygon.unsigned_area();

    let area_ratio = clipped_area / original_area;
    log_1(&format!("Area ratio: {}", area_ratio).into());
    JsValue::from(area_ratio)
}

// Generates random trees for all strata in stand with jittered grid sampling
pub fn generate_random_trees_into_buffer(
    stand_p: &Polygon,
    clipped_p: &Polygon,
    strata: &TreeStrata,
    stand_number: f64,
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
                    //if clipped_p.contains(&point!(x: pair[0], y: pair[1])) {
                    Tree::new(
                        stand_number,
                        stratum.tree_species,
                        stratum.mean_height,
                        (pair[0], pair[1], 0.0),
                    )
                    //}
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
                tree.stand_number(),
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

    // Create a shared buffer to store the generated trees
    //let buffer = SharedBuffer::new(max_tree_count as usize);

    // If there are stands in the bounding box, generate random trees for each stand
    if let Some(stands) = stands {
        let mut compartment_areas = Vec::new();
        let mut total_tree_count = 0;

        /*         let mut buffer_index = 0; */
        for stand in stands {
            let polygon = stand.computed_polygon.to_owned().unwrap();
            /*          let strata = stand.get_strata();
            let stand_number: f64 = stand.stand_basic_data.stand_number as f64; */

            // Clip the stand's polygon to the bounding box
            let intersected_polygons = polygon.intersection(bbox).0;
            let clipped_polygon = intersected_polygons
                .first()
                .expect("Intersection result should contain at least one polygon")
                .to_owned();

            // Generate trees and save them to the buffer if strata exist
            /* let mut tree_count = 0;
                       if let Some(strata) = strata {
                           tree_count = generate_random_trees_into_buffer(
                               &polygon,
                               &clipped_polygon,
                               &strata,
                               stand_number,
                               &buffer,
                               buffer_index,
                           );
                           buffer_index += tree_count;
                           log_1(
                               &format!(
                                   "Generated {} trees for stand {}",
                                   tree_count, stand.stand_basic_data.stand_number
                               )
                               .into(),
                           );
                       }
                       total_tree_count += tree_count;
            */

            total_tree_count += stand.summary_stem_count().unwrap_or(0) as usize;
            // Add to the compartment areas list
            compartment_areas.push(CompartmentArea {
                stand_number: stand.stand_basic_data.stand_number.to_string(),
                polygon: clipped_polygon,
            });
        }
        return Some((compartment_areas, max_tree_count as usize, total_tree_count));
        // Get a slice of the buffer
        /* let buffer_slice: &[f64] =
        unsafe { std::slice::from_raw_parts(buffer.ptr(), buffer.len()) }; */
    } else {
        None
    }
    // Log the buffer contents
    /*         log_1(&"Bounding box contains:".into());
    for (i, value) in buffer_slice.iter().enumerate() {
        if i % 6 == 0 && buffer_slice[i + 3] != 0.0 {
            let buffer_info = format!(
                "Tree {}: stand: {}, x = {}, y = {}, species = {}, height = {}, status = {}",
                i / 6,
                buffer_slice[i],
                buffer_slice[i + 1],
                buffer_slice[i + 2],
                buffer_slice[i + 3],
                buffer_slice[i + 4],
                buffer_slice[i + 5]
            );
            log_1(&buffer_info.into());
        }
    }
     */
    /*         let mut buffer_pointer_string = format!("{:p}", buffer.ptr());
    buffer_pointer_string = (&buffer_pointer_string[2..]).to_string(); */

    /*         log_1(&format!("Hexadecimal Buffer pointer in rust: {}", buffer_pointer_string).into());
    let buffer_pointer = hexadecimal_to_decimal(&buffer_pointer_string).unwrap();
    log_1(&format!("Decimal Buffer pointer in rust: {}", buffer_pointer).into()); */
    /*  (compartment_areas, max_tree_count as usize, total_tree_count)
    } else {
        (vec![], 0, 0, 0)
    } */
}

/* pub fn hexadecimal_to_decimal(hexadecimal_str: &str) -> Result<u64, &'static str> {
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
 */
#[wasm_bindgen]
pub fn empty_function(xml_content: String) {}
