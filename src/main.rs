use geo::BooleanOps;
use geo_points::forest_property::forest_property_data::ForestPropertyData;
use geo_points::main_functions::{
    create_geo_json_from_coords, 
    draw_selected_stand, 
    draw_stands_in_bbox, 
    get_bounding_box_of_map, 
    random_bbox, 
    save_geojson
};
use geo_points::geometry_utils::get_min_max_coordinates;
use geo_points::requests::{fetch_buildings, buildings_as_polygons, fetch_roads};
use geojson::{FeatureCollection, GeoJson};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>>{
    let mut bbox = get_bounding_box_of_map();
    bbox = random_bbox(&bbox);
    println!("Bounding box: {:#?}", bbox);
    let empty_geojson = GeoJson::FeatureCollection(FeatureCollection {
        bbox: None,
        features: vec![],
        foreign_members: None,
    });
    
    let mut buildings_geojson = empty_geojson.clone();
    let mut buildings = Vec::new();
    match fetch_buildings(&bbox) {
        Ok(geojson) => {
            buildings = buildings_as_polygons(&geojson)?;
            println!("Fetched {} buildings", buildings.len());

            // Exclude buildings from the bounding box
            for building in buildings.iter() {
                bbox = bbox.difference(building).0.first().unwrap().to_owned();
            }

            buildings_geojson = geojson;
        }
        Err(e) => {
            eprintln!("Failed to fetch buildings: {}", e);
            return Err(Box::new(e));
        }
    }

    let mut roads_geojson = empty_geojson.clone();
    match fetch_roads(&bbox) {
        Ok(geojson) => {
            roads_geojson = geojson;
        }
        Err(e) => {
            eprintln!("Failed to fetch roads: {}", e);
            return Err(Box::new(e));
        }
    }

    let(min_x, max_x, min_y, max_y) = get_min_max_coordinates(&bbox);
    let property = ForestPropertyData::from_xml_file("forestpropertydata.xml");

    match create_geo_json_from_coords(min_x, max_x, min_y, max_y, &property, &buildings_geojson, &roads_geojson) {
        Ok(geojson) => {
            let filename = "stands_in_bbox.geojson";
            save_geojson(&geojson, filename);
            println!("GeoJson data successfully saved to {}", filename);
        }
        Err(e) => {
            eprintln!("Failed to create GeoJson data: {}", e);
            return Err(e); 
        }
    }

    println!("------------------------------------------------------------");
    let map_image = draw_stands_in_bbox(&bbox, &property, &buildings);  
    map_image
        .img()
        .save("stands_in_bbox_image.png")
        .expect("Failed to save image");
    println!("Image saved as 'stands_in_bbox_image.png'");

    println!("------------------------------------------------------------");
    let selected_image = draw_selected_stand(&property);
    selected_image
        .img()
        .save("selected_stand_image.png")
        .expect("Failed to save image");
    println!("Image saved as {}", "selected_stand_image.png");

    Ok(())
}
