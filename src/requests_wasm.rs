use crate::forest_property;
use crate::forest_property::compartment::{get_compartment_areas_in_bounding_box, get_compartments_in_bounding_box};
use crate::forest_property::forest_property_data::{ForestPropertyData, RealEstate, ForestPropertyDataSchema};
use crate::forest_property::geometry::PolygonGeometry;
use crate::forest_property::tree::Tree;
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::forest_property::trees::Trees;
use crate::geojson_utils::{all_compartment_areas_to_geojson, geojson_to_polygons, get_geojson_from_url, roads_geojson_to_linestrings, water_geojson_to_polygons, FetchError};
use crate::geometry_utils::get_coords_of_map;

use geo::{coord, point, BooleanOps, Closest, Contains, LineString, Point, Polygon};
use geo::algorithm::closest_point::ClosestPoint;
use geo::algorithm::haversine_distance::HaversineDistance;
use geojson::GeoJson;
use reqwest_wasm::Client;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::JsValue;
use web_sys::console::log_1;

const METERS_IN_ONE_DEGREE_LAT: f64 = 111_320.0;

fn meters_to_degrees_lat(meters: f64) -> f64 {
    meters / METERS_IN_ONE_DEGREE_LAT
}

fn meters_to_degrees_lon(meters: f64, latitude: f64) -> f64 {
    let latitude_radians = latitude.to_radians();
    let meters_per_degree_lon = METERS_IN_ONE_DEGREE_LAT * latitude_radians.cos();
    meters / meters_per_degree_lon
}

const THRESHOLD: f64 = 5.0; // 5 meters in meters

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum OperationType {
    Thinning(f64),
    Cutting(f64),
    Simulation(TreeStrata),
}

impl OperationType {
    pub fn new_thinning(volume: f64) -> Self {
        OperationType::Thinning(volume)
    }

    pub fn new_cutting(volume: f64) -> Self {
        OperationType::Cutting(volume)
    }

    pub fn new_simulation(strata: TreeStrata) -> Self {
        OperationType::Simulation(strata)
    }

    pub fn check_simulation(&self) -> bool {
        match self {
            OperationType::Simulation(_tree_strata_vec) => {
                // The variant is `Simulation`, and it contains `TreeStrata`
                true
            }
            _ => {
                // It's not the `Simulation` variant
                false
            }
        }
    }
    
    pub fn get_simulation_strata(&self) -> TreeStrata {
        match self {
            OperationType::Simulation(tree_strata) => {
                // The variant is `Simulation`, and it contains `Vec<TreeStrata>`
                tree_strata.to_owned()
            }
            _ => {
                // It's not the `Simulation` variant
                TreeStrata::new(Vec::new())
            }
        }
    }
    
    pub fn get_cutting_volume(&self) -> f64 {
        match self {
            OperationType::Thinning(volume) => *volume,
            OperationType::Cutting(volume) => *volume,
            _ => 0.0,
        }
    }

    pub fn get_type(&self) -> String {
        match self {
            OperationType::Thinning(_) => "Thinning".to_string(),
            OperationType::Cutting(_) => "Cutting".to_string(),
            OperationType::Simulation(_) => "Simulation".to_string(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[wasm_bindgen]
struct Operation {
    operation_type: OperationType,
    cutting_areas: Vec<Polygon>,
}

impl Operation {
    pub fn new(
        operation_type: OperationType, 
        cutting_areas: Vec<Polygon>
    ) -> Self {

        Operation {
            operation_type,
            cutting_areas,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[wasm_bindgen]
struct StandOperation {
    id: u32,
    stand_id: u32,
    operation: Operation,
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

// Fetches ForestPropertyData from the given url
pub async fn get_xml_from_url(url: String) -> Result<ForestPropertyData, FetchError> {
    let client = Client::new();

    let response = client.get(&url).send().await.map_err(FetchError::from)?;
    let xml = response.text().await.map_err(FetchError::from)?;
    let property = ForestPropertyData::parse_from_str(xml.as_str());

    Ok(property)
}

#[wasm_bindgen]
impl VirtualForest {
    #[wasm_bindgen(constructor)]
    pub fn new(xml: &str) -> Self {
        let property = ForestPropertyData::from_xml_str(xml);
        //let property = ForestPropertyData::default();

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

    // Method to fetch forest property data from a polygon string and update the VirtualForest instance
    #[wasm_bindgen]
    pub async fn update_property_data_from_polygon_string(
        &mut self, 
        polygon_string: String
    ) -> Result<(), JsValue> {

        // Construct the URL based on the polygon string
        let url = format!(
            "https://avoin.metsakeskus.fi/rest/mvrest/FRStandData/v1/ByPolygon?wktPolygon=POLYGON%20(({}))&stdVersion=MV1.9",
            polygon_string
        );

        // Fetch the data and update the property field
        match get_xml_from_url(url).await {
            Ok(property) => {
                // Update the VirtualForest's property field with the fetched data
                self.property = property;
                Ok(())
            }
            Err(_) => {
                log_1(&format!("Failed to fetch the ForestPropertyData from the URL").into());
                Err(JsValue::from_str("Failed to fetch the ForestPropertyData"))
            }
        }
    }

    // Gets the infrastructure data (buildings and roads) in the bounding box
    pub async fn get_infrastructure(
        &mut self, 
        xml: &str
    ) {
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

                    let new_polygons = water_geojson_to_polygons(&geojson);
        
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

    fn get_url_water(
        body_type: &str, 
        west: f64, 
        south: f64, 
        east: f64, 
        north: f64
    ) -> String {

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
            .clone()
            .unwrap()
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
    pub fn set_realestate(
        &mut self, 
        index: u32
    ) {
        self.selected_realestate = index;
    }

    #[wasm_bindgen]
    pub fn get_selected_realestate(&self) -> Result<RealEstateValue, JsValue> {
        if let Some((index, result)) = self
            .property
            .real_estates
            .clone()
            .unwrap()
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
        if let Some((_index, result)) = self
            .property
            .real_estates
            .clone()
            .unwrap()
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
    pub fn get_stand_by_id(
        &self, 
        id: String
    ) -> Result<JsValue, JsValue> {

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

    fn is_point_within_threshold(
        road_point: &Point<f64>, 
        tree: &mut Tree, 
        threshold_in_meters: f64
    ) -> bool {

        let point = point!(x: tree.position().0, y: tree.position().1);

        // Calculate the minimum distance from the point to any segment of the line
        let min_distance = road_point.haversine_distance(&point);
    
        // Check if the minimum distance is within the threshold (e.g., 5 meters)
        min_distance <= threshold_in_meters
    }

    #[wasm_bindgen]
    pub fn generate_trees_bbox(
        &self, 
        min_x: f64, 
        max_x: f64, 
        min_y: f64, 
        max_y: f64
    ) -> Trees {

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
        let mut compartments = get_compartments_in_bounding_box(stands, &bbox);
 
        if let Some(so) = self.stand_operations.iter().find(|so| so.active_operation == 1) {
            log_1(&format!("Found active operation").into());
            let stand_id = so.stand_id.to_string();
            let op_type = so.operation.operation_type.clone();
            let areas = so.operation.cutting_areas.clone();
     
            if let Some(compartment) = compartments.iter_mut().find(|comp| comp.stand_id == stand_id) {
                // Operate on the compartment with the specified operation type and areas
                log_1(&format!("Operating on compartment with stand number {}", stand_id).into());
                compartment.operate_compartment(op_type, areas);
            }
        }

        let road_lines = self.roads.clone().unwrap_or(vec![]);
        let mut trees = Trees::new(Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        
        compartments.iter_mut().for_each(|compartment| {
            compartment.trees.iter_mut().for_each(|tree| {
                // Move the tree from the road if it is within the threshold
/*                 let (x, y, _) = tree.position();
                let tree_point = point!(x: x, y: y);

                if let Some((mut road_point, road_line)) = road_lines.iter().filter_map(|rl| {
                    // Find the closest point on the road line to the tree point
                    match rl.closest_point(&tree_point) {
                        // If the point is within the threshold, return the point and the road line
                        Closest::SinglePoint(pt) if Self::is_point_within_threshold(&pt, tree, THRESHOLD) => Some((pt, rl)),
                        _ => None,
                    }
                }).next() {
                    Self::move_point_from_road(tree, &mut road_point, road_line, &compartment.polygon);
                } */
                
                trees.insert(*tree);
            });
        });

        log_1(&format!("Done generating trees!").into());

        return trees;
    }

    fn move_point_from_road(
        tree: &mut Tree, 
        road_point: &mut Point<f64>, 
        road_line: &LineString<f64>, 
        compartment: &Polygon<f64>
    ) {
        let (mut x, mut y, _) = tree.position();
        
        let five_meters_x = meters_to_degrees_lon(5.0, y);
        let five_meters_y = meters_to_degrees_lat(5.0);
        
        let dx = x - road_point.x();
        let dist_x = dx.abs();

        // Avoid division by zero
        if dist_x < 1e-9 {
            for dx_multiplier in [-1.0, 1.0] {
                let new_x = road_point.x() + (dx_multiplier * five_meters_x);

                if compartment.contains(&point!(x: new_x, y: y)) {
                    x = new_x;

                    // Update the road's point after moving on x-axis
                    if let Closest::SinglePoint(pt) = road_line.closest_point(&point!(x: x, y: y)) {
                        *road_point = pt;
                    }

                    log_1(&format!("Moved tree on x-axis from road").into());
                    break;
                }
            }
        } else if dist_x < five_meters_x {
            let scale_x = five_meters_x / dist_x;
            x = road_point.x() + dx * scale_x;

            // Update the road's point after moving on x-axis
            if let Closest::SinglePoint(pt) = road_line.closest_point(&point!(x: x, y: y)) {
                *road_point = pt;
            }
        }

        let dy = y - road_point.y();
        let dist_y = dy.abs();

        if dist_y < 1e-9 {
            for dy_multiplier in [-1.0, 1.0] {
                let new_y = road_point.y() + (dy_multiplier * five_meters_y);

                if compartment.contains(&point!(x: x, y: new_y)) {
                    y = new_y; 

                    log_1(&format!("Moved tree on y-axis from road").into());
                    break;
                }
            }
        } else if dist_y < five_meters_y {
            let scale_y = five_meters_y / dist_y;
            y = road_point.y() + dy * scale_y;
        }

        tree.set_position((x, y, 0.0));
    }
    
    #[wasm_bindgen]
    pub fn set_operation(
        &mut self, 
        stand_id: u32, 
        operation_name: u32, 
        area_polygon: String, 
        cutting_volume: Option<f64>, 
        new_strata: JsValue
    ) {
        let polygon_geometry: PolygonGeometry = serde_json::from_str(&area_polygon)
            .map_err(|e| JsValue::from_str(&e.to_string())).expect("REASON");

        let polygon = polygon_geometry.polygon_property.polygon.clone();
        let polygon = polygon.to_geo_polygon();
 
        let mut polygons = Vec::new();
        polygons.push(polygon.clone());

        let operation_type = match operation_name {
            1 => OperationType::Cutting(cutting_volume.unwrap_or(0.0)),
            2 => OperationType::Thinning(cutting_volume.unwrap_or(0.0)),
            3 => {
                let strata = if new_strata.is_null() || new_strata.is_undefined() {
                    TreeStrata::new(Vec::new())
                } else {
                    let json_string = new_strata.as_string().unwrap(); // Get JSON string from JsValue
                    serde_json::from_str(&json_string).expect("REASON") // Deserialize into TreeStrata
                };
                OperationType::Simulation(strata)
            },
            _ => OperationType::Cutting(0.0),
        };

        let operation = Operation::new(operation_type, polygons);
        let stand_operation = StandOperation {
            id: self.stand_operations.len() as u32 + 1,
            stand_id,
            operation,
            active_operation: 1,
        };
    
        self.stand_operations.push(stand_operation);
    }
    
    
    // Fetches GeoJSON data from the given bounding box and XML content
    // Returns a GeoJson of stand polygons, building polygons, and a roads multi-linestring
    // Buildings, water bodies, and roads have a property "type" with values "building", "water", and "roads" respectively
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

        // Exclude water bodies from the bounding box
        for water_body in self.water.clone().unwrap().iter() {
            bbox = bbox.difference(water_body).0.first().unwrap().to_owned();
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