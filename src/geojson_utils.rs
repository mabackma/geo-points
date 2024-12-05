use crate::{forest_property::{compartment::{Compartment, CompartmentArea}, tree::Tree}, geometry_utils::get_min_max_coordinates};

use geo::{coord, LineString, Polygon};
use geojson::{Feature, FeatureCollection, GeoJson, Geometry as GeoJsonGeometry, JsonObject, JsonValue, Value};
use geojson::Error as GeoJsonError;
use reqwest::Error as ReqwestError;
use reqwest_wasm::Error as ReqwestWasmError;
use serde_json::Error as SerdeJsonError;
use std::fmt;
use reqwest_wasm::Client;

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

// Function to convert a Polygon into a GeoJSON Feature
fn convert_polygon_to_feature(polygon: &Polygon<f64>, property: Option<&str>) -> Feature {
    let exterior_coords: Vec<Vec<f64>> = polygon.exterior().points()
        .map(|point| vec![point.x(), point.y()])
        .collect();

    let geometry = GeoJsonGeometry {
        bbox: None,
        value: Value::Polygon(vec![exterior_coords]),
        foreign_members: None,
    };

    if let Some(property) = property {
        let mut properties = JsonObject::new();
        properties.insert("type".to_string(), JsonValue::from(property));
        Feature {
            geometry: Some(geometry),
            properties: Some(properties),
            id: None,
            bbox: None,
            foreign_members: None,
        }
    } else {
        Feature {
            geometry: Some(geometry),
            properties: None,
            id: None,
            bbox: None,
            foreign_members: None,
        }
    }
}

fn convert_linestrings_to_feature(line_strings: &Vec<LineString>) -> Feature {
    let mut line_strings_coords = Vec::new();
    for line_string in line_strings {
        let coords: Vec<Vec<f64>> = line_string.points()
            .map(|point| vec![point.x(), point.y()])
            .collect();
        line_strings_coords.push(coords);
    }

    let geometry = GeoJsonGeometry {
        bbox: None,
        value: Value::MultiLineString(line_strings_coords),
        foreign_members: None,
    };

    let mut properties = JsonObject::new();
    properties.insert("type".to_string(), JsonValue::from("roads".to_string()));

    Feature {
        geometry: Some(geometry),
        properties: Some(properties),
        id: None,
        bbox: None,
        foreign_members: None,
    }
}

// Function to convert a Tree into a GeoJSON Feature
pub fn convert_tree_to_feature(tree: &Tree) -> Feature {
    let point = vec![tree.position().0, tree.position().1];
    let point_geometry = GeoJsonGeometry {
        bbox: None,
        value: Value::Point(point),
        foreign_members: None,
    };

    let mut properties = serde_json::Map::new();
    properties.insert("species".to_string(), serde_json::json!(tree.species()));
    properties.insert("stand_id".to_string(), serde_json::json!(tree.stand_id()));
    properties.insert("status".to_string(), serde_json::json!(tree.tree_status()));

    Feature {
        geometry: Some(point_geometry),
        properties: Some(properties),
        id: None,
        bbox: None,
        foreign_members: None,
    }
}

pub fn all_compartments_to_geojson(
        compartments: Vec<Compartment>,
        buildings: &GeoJson, 
        roads: &GeoJson) -> GeoJson {
        
    let mut all_features = Vec::new();

    for compartment in compartments {        
        // Get the trees within the clipped polygon
        let (min_x, max_x, min_y, max_y) = get_min_max_coordinates(&compartment.polygon);
        let trees = compartment.trees_in_bounding_box(min_x, max_x, min_y, max_y);

        // Convert the compartment (polygon) to a GeoJSON feature
        let polygon_feature = convert_polygon_to_feature(&compartment.polygon, None);
        let tree_features: Vec<Feature> = trees.iter().map(|tree| convert_tree_to_feature(tree)).collect();

        // Add the polygon feature and tree features to the list
        all_features.push(polygon_feature);
        all_features.extend(tree_features);
    }

    // Add building features to the list, ensuring the GeoJson is a FeatureCollection
    if let GeoJson::FeatureCollection(building_collection) = buildings {
        println!("Added buildings to geojson: {}", building_collection.features.len());
        for building_feature in &building_collection.features {
            all_features.push(building_feature.clone());
        }
    } else {
        println!("Buildings GeoJson is not a FeatureCollection");
    }

    // Add road features to the list, ensuring the GeoJson is a FeatureCollection
    if let GeoJson::FeatureCollection(road_collection) = roads {
        println!("Added roads to geojson: {}", road_collection.features.len());
        for road_feature in &road_collection.features {
            all_features.push(road_feature.clone());
        }
    } else {
        println!("Roads GeoJson is not a FeatureCollection");
    }

    // Create the GeoJSON FeatureCollection
    let feature_collection = FeatureCollection {
        features: all_features,
        bbox: None,
        foreign_members: None,
    };

    // Create a GeoJson object
    let geojson = GeoJson::FeatureCollection(feature_collection);

    geojson
}

pub fn all_compartment_areas_to_geojson(
    compartment_areas: Vec<CompartmentArea>,
    buildings: &Vec<Polygon>, 
    roads: &Vec<LineString>,
    water: &Vec<Polygon>) -> GeoJson {
    
    let mut all_features = Vec::new();

    for compartment_area in compartment_areas {
        let polygon_feature = convert_polygon_to_feature(&compartment_area.polygon, None);
        all_features.push(polygon_feature);
    }

    for building in buildings.iter() {
        let building_feature = convert_polygon_to_feature(building, Some("building"));
        all_features.push(building_feature);
    }

    for water_body in water.iter() {
        let water_feature = convert_polygon_to_feature(water_body, Some("water"));
        all_features.push(water_feature);
    }

    let road_feature = convert_linestrings_to_feature(&roads);
    all_features.push(road_feature);


    // Create the GeoJSON FeatureCollection
    let feature_collection = FeatureCollection {
        features: all_features,
        bbox: None,
        foreign_members: None,
    };

    // Create a GeoJson object
    let geojson = GeoJson::FeatureCollection(feature_collection);

    geojson
}

pub fn polygon_to_geojson(polygon: &Polygon<f64>, trees: &Vec<Tree>) -> GeoJson {
    let mut all_features = Vec::new();

    // Convert the compartment (polygon) to a GeoJSON feature
    let polygon_feature = convert_polygon_to_feature(&polygon, None);
    let tree_features: Vec<Feature> = trees.iter().map(|tree| convert_tree_to_feature(tree)).collect();

    // Add the polygon feature and tree features to the list
    all_features.push(polygon_feature);
    all_features.extend(tree_features);

    // Create the GeoJSON FeatureCollection
    let feature_collection = FeatureCollection {
        features: all_features,
        bbox: None,
        foreign_members: None,
    };

    // Return a GeoJson object
    GeoJson::FeatureCollection(feature_collection)
}

// Fetches GeoJSON data from the given url
pub async fn get_geojson_from_url(url: String) -> Result<GeoJson, FetchError> {
    let client = Client::new();

    let response = client.get(&url).send().await.map_err(FetchError::from)?;
    let text = response.text().await.map_err(FetchError::from)?;
    let geojson: GeoJson = serde_json::from_str(&text).map_err(FetchError::from)?;

    Ok(geojson)
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

pub fn water_geojson_to_polygons(geojson: &GeoJson) -> Vec<Polygon<f64>> {
    let mut polygons = Vec::new();

    if let GeoJson::FeatureCollection(collection) = geojson {
        for feature in collection.features.iter() {
            if let Some(geometry) = &feature.geometry {
                match &geometry.value {
                    Value::Polygon(polygon) => {
                        let exterior = polygon[0]
                            .iter()
                            .filter_map(|point| {
                                if point[2] >= 0.0 {
                                    Some(coord!(x: point[0], y: point[1]))
                                } else {
                                    None
                                }
                            })
                            .collect::<Vec<_>>();
                        let poly = Polygon::new(LineString::from(exterior), vec![]);
                        polygons.push(poly);
                    }
                    _ => {
                        continue;
                    }
                }
            }
        }
    }

    polygons
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