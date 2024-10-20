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

use geo::{coord, point, Area, BooleanOps, Contains, EuclideanDistance, Line, LineString, Polygon};
use geojson::{Error as GeoJsonError, JsonObject};
use geojson::{GeoJson, Value};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use reqwest::Error as ReqwestError;
use reqwest_wasm::Client;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value as SerdeJsonValue};
use serde_wasm_bindgen;
use web_sys::js_sys::Math::abs;
use std::error::Error;
use std::fmt;
use std::slice::from_raw_parts;
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::{JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
use web_sys::console::log_1;
use web_sys::js_sys::{Float32Array, JsString};
use reqwest_wasm::Error as ReqwestWasmError;
use serde_json::Error as SerdeJsonError;

#[derive(Debug)]
pub enum FetchError {
    Reqwest(ReqwestError),
    GeoJson(GeoJsonError),
    ReqwestWasm(ReqwestWasmError),
    SerdeJson(SerdeJsonError),
}

impl fmt::Display for FetchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FetchError::Reqwest(err) => write!(f, "Reqwest error: {}", err),
            FetchError::GeoJson(err) => write!(f, "GeoJson error: {}", err),
            FetchError::ReqwestWasm(err) => write!(f, "Reqwest WASM error: {}", err),
            FetchError::SerdeJson(err) => write!(f, "Serde JSON error: {}", err),

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

impl From<ReqwestWasmError> for FetchError {
    fn from(err: ReqwestWasmError) -> Self {
        FetchError::ReqwestWasm(err)
    }
}

impl From<SerdeJsonError> for FetchError { 
    fn from(err: SerdeJsonError) -> Self {
        FetchError::SerdeJson(err)
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

#[derive(Serialize, Deserialize, Clone, PartialEq)]
enum OperationType {
    Thinning,
    Cutting,
    Simulation,
}

#[derive(Serialize, Deserialize, Clone, PartialEq)]
#[wasm_bindgen]
struct Operation {
    operation_type: OperationType,
    cutting_areas: Vec<Polygon>,
}

#[derive(Serialize, Deserialize, Clone, PartialEq)]
#[wasm_bindgen]
struct StandOperation {
    id: u32,
    stand_id: u32,
    operations: Vec<TreeStrata>,
    active_operation: u32,
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone, PartialEq)]
pub struct VirtualForest {
    property: ForestPropertyData,
    selected_realestate: u32,
    stand_operations: Vec<StandOperation>,
    retention_zones: Vec<Polygon>,
    roads: Option<Vec<LineString>>,
    buildings: Option<Vec<Polygon>>,
    water: Option<Vec<Polygon>>,
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
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
            buildings: None,
        }
    }

    // Gets the infrastructure data (buildings and roads) in the bounding box
    pub async fn get_infrastructure(&mut self, xml: &str) {
        let (min_x, max_x, min_y, max_y) = get_coords_of_map(xml);

        let west = min_x;
        let south = min_y;
        let east = max_x;
        let north = max_y;
    
        // Get buildings in the bounding box
        let url_buildings = format!(
            "https://metne-test.onrender.com/geoserver/mml/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=mml:rakennus&maxFeatures=2000&outputFormat=application%2Fjson&BBOX={},{},{},{},EPSG:4326&srsName=EPSG:4326",
            west, south, east, north
        );

        match get_geojson_from_url(url_buildings).await {
            Ok(geojson) => {
                self.buildings = Some(geojson_to_polygons(&geojson));
            }
            Err(_) => {
                self.buildings = None;
            }
        }

        // Get roads in the bounding box
        let url_roads = format!(
            "https://metne-test.onrender.com/geoserver/mml/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=mml:tieviiva&bbox={},{},{},{},EPSG:4326&srsName=EPSG:4326&outputFormat=application/json",
            west, south, east, north
        );

        match get_geojson_from_url(url_roads).await {
            Ok(geojson) => {
                self.roads = Some(roads_geojson_to_linestrings(&geojson));
            }
            Err(_) => {
                self.roads = None;
            }
        }
        
        // Get water bodies in the bounding box
        for body_type in ["allas", "jarvi", "meri", "virtavesialue"].iter() {
            match get_geojson_from_url(Self::get_url_water(body_type, west, south, east, north)).await {
                Ok(geojson) => {
                    let water_body_count = match &geojson {
                        GeoJson::FeatureCollection(collection) => collection.features.len(),
                        _ => 0,
                    };
                    log_1(&format!("Fetched water body of type: {} {}", body_type, water_body_count).into());

                    let new_polygons = geojson_to_polygons(&geojson);
        
                    // Append the new bodies of water to the existing ones
                    if let Some(existing_polygons) = &mut self.water {
                        existing_polygons.extend(new_polygons);
                    } else {
                        self.water = Some(new_polygons);
                    }
                }
                Err(_) => {
                    log_1(&format!("Failed to fetch water body of type: {}", body_type).into());
                }
            }
        }
    }

    fn get_url_water(body_type: &str, west: f64, south: f64, east: f64, north: f64) -> String {
        format!(
            "https://metne-test.onrender.com/geoserver/mml/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=mml:{}&bbox={},{},{},{},EPSG:4326&srsName=EPSG:4326&outputFormat=application/json",
            body_type, west, south, east, north
        )
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
        let property: Self = serde_json::from_str(json).expect("Error parsing json file");

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

                     /*    let l = LineString::new(vec![coord! {x: 10.0,y:10.0},coord! {x: 20.0,y:10.0}]);

                        let remove = l.points().find(|p| {

                            abs( p.x() - x ) > 0.0005 || abs(p.y() - y) > 0.0005

                        });


                        let p = point!(x: x,y: y);

                        let distance = p.euclidean_distance(&l); */

                        // let isOk = l.euclidean_distance(&point!(x: x as usize,y: y as usize)) > 5.0 as usize;



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
    // Returns a GeoJson of stand polygons, building polygons, and a roads multi-linestring
    // Buildings and roads have a property "type" with values "building" and "roads" respectively
    #[wasm_bindgen]
    pub async fn geo_json_from_coords(
        & mut self,
        min_x: f64,
        max_x: f64,
        min_y: f64,
        max_y: f64,
    ) -> Result<JsValue, JsValue> {

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

        let buildings_count = self.buildings.clone().unwrap().len();
        log_1(&format!("{} buildings in the bounding box", buildings_count).into());
        log_1(&format!("{} roads in the bounding box", self.roads.clone().unwrap().len()).into());
        log_1(&format!("{} water bodies in the bounding box", self.water.clone().unwrap().len()).into());

        // Exclude buildings from the bounding box
        for building in self.buildings.clone().unwrap().iter() {
            bbox = bbox.difference(building).0.first().unwrap().to_owned();
        }

        // Get the ForestPropertyData and stands
        let stands = self._get_selected_realestate().unwrap().get_stands();

        // Get compartment areas in the bounding box and convert them to GeoJSON
        if let Some(compartment_areas) = get_compartment_areas_in_bounding_box(stands, &bbox) {

            let geojson = all_compartment_areas_to_geojson(
                compartment_areas.0,
                &self.buildings.clone().unwrap_or(vec![]),
                &self.roads.clone().unwrap_or(vec![]),
                &self.water.clone().unwrap_or(vec![]),
            );

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
                    Tree::new(
                        stand_number,
                        stratum.tree_species,
                        stratum.mean_height,
                        (pair[0], pair[1], 0.0),
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

// Get the bounding box of the whole map
pub fn get_coords_of_map(xml: &str) -> (f64, f64, f64, f64) {
    let property = ForestPropertyData::from_xml_str(xml);
    let mut all_stands = property.real_estates.real_estate[0].get_stands();

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

// Fetches GeoJSON data from the given url
pub async fn get_geojson_from_url(url: String) -> Result<GeoJson, FetchError> {
    let client = Client::new();

    let response = client.get(&url).send().await.map_err(FetchError::from)?;
    let text = response.text().await.map_err(FetchError::from)?;
    let geojson: GeoJson = serde_json::from_str(&text).map_err(FetchError::from)?;

    Ok(geojson)
}

pub fn roads_geojson_to_linestrings(geojson: &GeoJson) -> Vec<LineString<f64>> {
    let mut linestrings = Vec::new();

    if let GeoJson::FeatureCollection(collection) = geojson {
        for feature in collection.features.iter() {
            if let Some(geometry) = &feature.geometry {
                match &geometry.value {
                    Value::LineString(line_string) => {
                        let line = line_string
                            .iter()
                            .filter_map(|point| {
                                if point[2] >= 0.0 {
                                    Some(coord!(x: point[0], y: point[1]))
                                } else {
                                    None
                                }
                            })
                            .collect::<Vec<_>>();
                        linestrings.push(LineString::from(line));
                    }
                    _ => {
                        continue;
                    }
                }
            }
        }
    }

    linestrings
}

