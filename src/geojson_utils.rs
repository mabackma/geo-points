use crate::{forest_property::{compartment::{Compartment, CompartmentArea}, tree::Tree}, geometry_utils::get_min_max_coordinates};

use geo::Polygon;
use geojson::{Feature, FeatureCollection, GeoJson, Geometry as GeoJsonGeometry, Value};

// Function to convert a Polygon into a GeoJSON Feature
fn convert_polygon_to_feature(polygon: &Polygon<f64>) -> Feature {
    let exterior_coords: Vec<Vec<f64>> = polygon.exterior().points()
        .map(|point| vec![point.x(), point.y()])
        .collect();

    let geometry = GeoJsonGeometry {
        bbox: None,
        value: Value::Polygon(vec![exterior_coords]),
        foreign_members: None,
    };

    Feature {
        geometry: Some(geometry),
        properties: None,
        id: None,
        bbox: None,
        foreign_members: None,
    }
}

// Function to convert a Tree into a GeoJSON Feature
fn convert_tree_to_feature(tree: &Tree) -> Feature {
    let point = vec![tree.position().0, tree.position().1];
    let point_geometry = GeoJsonGeometry {
        bbox: None,
        value: Value::Point(point),
        foreign_members: None,
    };

    let mut properties = serde_json::Map::new();
    properties.insert("species".to_string(), serde_json::json!(tree.species()));

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
        let polygon_feature = convert_polygon_to_feature(&compartment.polygon);
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
    buildings: &GeoJson, 
    roads: &GeoJson) -> GeoJson {
    
    let mut all_features = Vec::new();

    for compartment_area in compartment_areas {
        // Convert the compartment (polygon) to a GeoJSON feature
        let polygon_feature = convert_polygon_to_feature(&compartment_area.polygon);

        // Add the polygon feature to the list
        all_features.push(polygon_feature);
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

pub fn polygon_to_geojson(polygon: &Polygon<f64>, trees: &Vec<Tree>) -> GeoJson {
    let mut all_features = Vec::new();

    // Convert the compartment (polygon) to a GeoJSON feature
    let polygon_feature = convert_polygon_to_feature(&polygon);
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