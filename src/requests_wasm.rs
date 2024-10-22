use crate::forest_property::compartment::{get_compartment_areas_in_bounding_box, get_compartments_in_bounding_box};
use crate::forest_property::forest_property_data::{ForestPropertyData, RealEstate};
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::geojson_utils::{all_compartment_areas_to_geojson, geojson_to_polygons, get_geojson_from_url, roads_geojson_to_linestrings, water_geojson_to_polygons};
use crate::geometry_utils::get_coords_of_map;

use geo::{coord, point, BooleanOps, Closest, Contains, EuclideanDistance, LineString, Point, Polygon};
use geo::algorithm::closest_point::ClosestPoint;
use geo::algorithm::haversine_distance::HaversineDistance;
use geojson::GeoJson;
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

    fn is_point_within_threshold(road_point: &Point<f64>, point: &Point<f64>, threshold_in_meters: f64) -> bool {
        // Calculate the minimum distance from the point to any segment of the line
        let min_distance = road_point.haversine_distance(point);
    
        // Check if the minimum distance is within the threshold (e.g., 5 meters)
        min_distance <= threshold_in_meters
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
        let road_lines = self.roads.clone().unwrap_or(vec![]);

        let trees = compartments
            .iter()
            .flat_map(|compartment| {
                compartment.trees
                    .iter()
                    .filter_map(|tree| {
                        let (mut x, mut y, z) = tree.position();
                        let specie = tree.species();
                        let height = tree.tree_height();
                        let status = tree.tree_status();
                        let stand_number = tree.stand_number();
                        
                        let tree_point = point!(x: x, y: y);

                        if let Some(road_point) = road_lines.iter().filter_map(|rl| {
                            // Find the closest point on the road line to the tree point
                            match rl.closest_point(&tree_point) {
                                // If the point is within the threshold, return the point
                                Closest::SinglePoint(pt) if Self::is_point_within_threshold(&pt, &tree_point, THRESHOLD) => Some(pt),
                                _ => None,
                            }
                        }).next() {
                            log_1(&format!("Moving point from road").into());
                            Self::move_point_from_road(&mut x, &mut y, &road_point, &compartment.polygon);
                        }

                        vec![
                            x as f32,
                            y as f32,
                            z as f32,
                            specie as f32,
                            height,
                            status as f32,
                            stand_number as f32,
                        ].into()  
                    })
            })
            .flatten()
            .collect::<Vec<f32>>();

        return trees;
    }

    fn move_point_from_road(x: &mut f64, y: &mut f64, road_point: &Point<f64>, compartment: &Polygon<f64>) {
        let (road_x, road_y) = road_point.x_y();
        let dx = *x - road_x;
        let dy = *y - road_y;
        let five_meters_x = meters_to_degrees_lon(5.0, *y);
        let five_meters_y = meters_to_degrees_lat(5.0);
        
        let dist = road_point.euclidean_distance(&point!(x: *x, y: *y));

        // Avoid division by zero
        if dist < 1e-9 {
            let movements = [
                (1.0, 1.0),  // Move +X +Y
                (1.0, -1.0), // Move +X -Y
                (-1.0, 1.0), // Move -X +Y
                (-1.0, -1.0) // Move -X -Y
            ];

            for (dx_multiplier, dy_multiplier) in &movements {
                let new_x = road_x + (dx_multiplier * five_meters_x);
                let new_y = road_y + (dy_multiplier * five_meters_y);

                // Check if the new position is within the compartment
                if compartment.contains(&point!(x: new_x, y: new_y)) {
                    *x = new_x; 
                    *y = new_y;
                    log_1(&format!("Moved tree from road").into());
                    break;
                }
            }
        } else {
            let scale_x = five_meters_x / dist;
            let scale_y = five_meters_y / dist;

            // Move the point 5 meters away from the road
            *x = road_x + dx * scale_x;
            *y = road_y + dy * scale_y;
        }
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
