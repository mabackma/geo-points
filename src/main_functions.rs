use std::fs::File;
use crate::geometry_utils::{generate_random_trees, get_min_max_coordinates};
use crate::geojson_utils::{polygon_to_geojson, all_compartments_to_geojson};
use crate::forest_property::compartment::get_compartments_in_bounding_box;
use crate::forest_property::forest_property_data::ForestPropertyData;
use crate::forest_property::image_processor::ImageProcessor;
use geo::{coord, Coord, LineString, Polygon};
use geojson::GeoJson;
use image::Rgb;
use rand::Rng;
use std::io::Write;
use std::time::Instant;
use std::error::Error;


// Get the bounding box of the whole map
pub fn get_bounding_box_of_map() -> Polygon<f64> {
    let property = ForestPropertyData::from_xml_file("forestpropertydata.xml");
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
    
    let bbox = geo::Polygon::new(
        LineString(vec![
            Coord { x: min_x, y: min_y },
            Coord { x: max_x, y: min_y },
            Coord { x: max_x, y: max_y },
            Coord { x: min_x, y: max_y },
            Coord { x: min_x, y: min_y },
        ]),
        vec![],
    );

    bbox
}

pub fn random_bbox(map_bbox: &Polygon<f64>) -> Polygon<f64> {
    let (min_x, max_x, min_y, max_y) = get_min_max_coordinates(&map_bbox);

    let mut rng = rand::thread_rng();

    let x1 = rng.gen_range(min_x..max_x);
    let y1 = rng.gen_range(min_y..max_y);
    let x2 = rng.gen_range(min_x..max_x);
    let y2 = rng.gen_range(min_y..max_y);

    let random_bbox = geo::Polygon::new(
        LineString(vec![
            Coord { x: x1, y: y1 },
            Coord { x: x2, y: y1 },
            Coord { x: x2, y: y2 },
            Coord { x: x1, y: y2 },
            Coord { x: x1, y: y1 },
        ]),
        vec![],
    );

    random_bbox
}

// Get color based on species number
fn get_color_by_species(number: u8) -> Rgb<u8> {
    match number {
        // Coniferous Trees (Shades of Orange and Red)
        1 => Rgb([255, 165, 0]),    // Orange - Mänty
        2 => Rgb([255, 0, 0]),      // Red - Kuusi
        10 => Rgb([255, 140, 0]),   // DarkOrange - Douglaskuusi
        11 => Rgb([255, 99, 71]),   // Tomato - Kataja
        12 => Rgb([255, 127, 80]),  // Coral - Kontortamänty
        16 => Rgb([178, 34, 34]),   // Firebrick - Mustakuusi
        19 => Rgb([205, 92, 92]),   // IndianRed - Pihta
        22 => Rgb([139, 0, 0]),     // DarkRed - Sembramänty
        23 => Rgb([233, 150, 122]), // DarkSalmon - Serbiankuusi
        30 => Rgb([250, 128, 114]), // Salmon - Havupuu

        // Deciduous Trees (Shades of Green and Blue)
        3 => Rgb([50, 205, 50]),    // LimeGreen - Rauduskoivu
        4 => Rgb([34, 139, 34]),    // ForestGreen - Hieskoivu
        5 => Rgb([107, 142, 35]),   // OliveDrab - Haapa
        6 => Rgb([143, 188, 143]),  // DarkSeaGreen - Harmaaleppä
        7 => Rgb([46, 139, 87]),    // SeaGreen - Tervaleppä
        9 => Rgb([32, 178, 170]),   // LightSeaGreen - Muu lehtipuu
        13 => Rgb([0, 128, 128]),   // Teal - Kynäjalava
        14 => Rgb([102, 205, 170]), // MediumAquamarine - Lehtikuusi
        15 => Rgb([60, 179, 113]),  // MediumSeaGreen - Metsälehmus
        17 => Rgb([152, 251, 152]), // PaleGreen - Paju
        18 => Rgb([0, 255, 127]),   // SpringGreen - Pihlaja
        20 => Rgb([0, 250, 154]),   // MediumSpringGreen - Raita
        21 => Rgb([144, 238, 144]), // LightGreen - Saarni
        24 => Rgb([85, 107, 47]),   // DarkOliveGreen - Tammi
        25 => Rgb([154, 205, 50]),  // YellowGreen - Tuomi
        26 => Rgb([0, 255, 0]),     // Lime - Vaahtera
        27 => Rgb([173, 216, 230]), // LightBlue - Visakoivu
        28 => Rgb([72, 209, 204]),  // MediumTurquoise - Vuorijalava
        29 => Rgb([64, 224, 208]),  // Turquoise - Lehtipuu

        // Default case for any unknown tree number
        _ => Rgb([0, 0, 0]), // Black for Unknown
    }
}

/* CREATES GEOJSON FROM COORDINATES OF BOUNDING BOX */
pub fn create_geo_json_from_coords(min_x: f64, max_x: f64, min_y: f64, max_y: f64, property: &ForestPropertyData, buildings_geojson: &GeoJson, roads_geojson: &GeoJson) -> Result<GeoJson, Box<dyn Error>>  {
    let start = Instant::now();

    let bbox = geo::Polygon::new(
        LineString(vec![
            coord!(x: min_x, y: min_y),
            coord!(x: max_x, y: min_y),
            coord!(x: max_x, y: max_y),
            coord!(x: min_x, y: max_y),
            coord!(x: min_x, y: min_y),
        ]),
        vec![],
    );

    let real_estate = property.real_estates.real_estate[0].clone();
    let stands = real_estate.get_stands();
    println!("Total stands: {:?}", stands.len());

    // Create compartments in the bounding box
    let compartments = get_compartments_in_bounding_box(stands, &bbox);
    println!("\nCompartments in bounding box: {:?}", compartments.len());

    let geojson = all_compartments_to_geojson(compartments, &buildings_geojson, &roads_geojson);

    let duration = start.elapsed();
    println!("\nTime elapsed in create_geo_json_for_bbox is: {:?}\n", duration);

    // Return all compartments and trees as GeoJson
    Ok(geojson)
}

pub fn draw_stands_in_bbox(bbox: &Polygon<f64>, property: &ForestPropertyData, buildings: &Vec<Polygon>) -> ImageProcessor {
    let start = Instant::now();

    let real_estate = property.real_estates.real_estate[0].clone();
    let stands = real_estate.get_stands();
    println!("Total stands: {:?}\n", stands.len());

    // Find compartments in the bounding box
    let compartments = get_compartments_in_bounding_box(stands, &bbox);
    println!("\nCompartments in bounding box: {:?}", compartments.len());

    let (min_x, max_x, min_y, max_y) = get_min_max_coordinates(&bbox);

    // Create an image processor with the desired image dimensions
    let img_width = ((max_x - min_x) * 100000.0) as u32;
    let img_height = ((max_y - min_y) * 100000.0) as u32;
    let mut image = ImageProcessor::new(img_width, img_height);

    let scale = ImageProcessor::create_scale(min_x, max_x, min_y, max_y, img_width, img_height);

    for compartment in compartments {
        let polygon = match compartment.clip_polygon_to_bounding_box(&bbox) {
            Some(polygon) => polygon,
            None => continue,
        };
        
        // Get the trees within the clipped polygon
        let (min_x, max_x, min_y, max_y) = get_min_max_coordinates(&polygon);
        let trees = compartment.trees_in_bounding_box(min_x, max_x, min_y, max_y);

        // Draw the clipped polygon
        let mapped_coordinates = image.map_coordinates_to_image(&polygon, &scale);
        image.draw_polygon_image(&mapped_coordinates, Rgb([0, 0, 255]));

        // Draw the trees
        for tree in trees {
            let point = coord! {x: tree.position().0, y: tree.position().1};
            let color = get_color_by_species(tree.species());
            image.draw_random_point(&scale, img_width, img_height, point, color);
        }
    }

    // Draw the buildings
    for building in buildings.iter() {
        let mapped_building = image.map_coordinates_to_image(&building, &scale);
        image.draw_polygon_image(&mapped_building, Rgb([255, 255, 255]));
    }

    let duration = start.elapsed();
    println!("\nTime elapsed in draw_stands_in_bbox is: {:?}", duration);

    image
}

/* ASKS USER FOR STAND AND DRAWS STAND. SAVES STAND TO GEOJSON */
pub fn draw_selected_stand(property: &ForestPropertyData) -> ImageProcessor {
    let mut stand = property.get_stand_cli();
    let polygon = stand.create_polygon();

    // Create an image for the polygon and random points
    let img_width = 800;
    let img_height = 600;
    let mut image = ImageProcessor::new(img_width, img_height);

    // Get the minimum and maximum x and y coordinates of the polygon
    let (min_x, max_x, min_y, max_y) = get_min_max_coordinates(&polygon);
    let scale = ImageProcessor::create_scale(min_x, max_x, min_y, max_y, img_width, img_height);

    // Map polygon coordinates to image
    let mapped_coordinates = image.map_coordinates_to_image(&polygon, &scale);
    image.draw_polygon_image(&mapped_coordinates, Rgb([0, 0, 255]));

    let summary_stem_count = stand.summary_stem_count();
    let strata = stand.get_strata().expect("No treeStrata/stratums found");
    let random_trees = generate_random_trees(&polygon, &strata, 1.0);

    // Convert the Polygon and the trees to GeoJSON
    let geojson = polygon_to_geojson(&polygon, &random_trees);

    // Serialize GeoJson to a String
    let geojson_string = serde_json::to_string_pretty(&geojson).expect("Failed to serialize GeoJson");
    
    // Write GeoJson to a file
    let mut file = File::create("selected_stand.geojson").expect("Failed to create file");
    file.write_all(geojson_string.as_bytes()).expect("Failed to write to file");
    
    if stand.stem_count_in_stratum() {
        println!("\nStem count is in individual stratum");

        // Draw the random points
        for tree in random_trees {
            let point = coord! {x: tree.position().0, y: tree.position().1};
            let color = get_color_by_species(tree.species());
            image.draw_random_point(&scale, img_width, img_height, point, color);
        }
    } else {
        println!("Stem count is not in any individual stratum. Drawing random points based on tree stand summary.");

        // Draw the random points
        for tree in random_trees {
            let point = coord! {x: tree.position().0, y: tree.position().1};
            image.draw_random_point(&scale, img_width, img_height, point, Rgb([255, 0, 0])) // Draw points in red
        }
    }
    println!("\nTotal stem count: {:?}", summary_stem_count);

    image
}

// Function to save a GeoJson object to a file
pub fn save_geojson(geojson: &GeoJson, filename: &str) {
    // Serialize the GeoJson object to a string
    let geojson_string = serde_json::to_string_pretty(&geojson).expect("Failed to serialize GeoJson");

    // Save the GeoJSON string to a file
    let mut file = File::create(filename).expect("Failed to create file");
    file.write_all(geojson_string.as_bytes()).expect("Failed to write to file");

    println!("GeoJSON saved to {}", filename);
}
