use geo::{Coord, LineString, Polygon};
use serde::{Deserialize, Serialize};
use crate::forest_property::tree_stand_data::TreeStrata;
use crate::forest_property::forest_property_data::{ TreeStandDataDate, TreeStratum};
use crate::forest_property::forest_property_data::{Operations, SpecialFeatures, StandBasicData, TreeStandData};
use crate::projection::{Projection, CRS};


#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Stands {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "Stand")]
    pub stand: Vec<Stand>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Stand {
    #[serde(rename = "@id")]
    pub id: String,
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "StandBasicData")]
    pub stand_basic_data: StandBasicData,
    #[serde(rename = "SpecialFeatures")]
    pub special_features: Option<SpecialFeatures>,
    #[serde(rename = "Operations")]
    pub operations: Option<Operations>,
    #[serde(rename = "TreeStandData")]
    pub tree_stand_data: Option<TreeStandData>,
    #[serde(skip_serializing, skip_deserializing)]
    pub computed_polygon: Option<Polygon>,
    #[serde(skip_serializing, skip_deserializing)]
    pub proj: Projection,
}

impl Stand {
    pub fn compute_polygon(&mut self) -> &Self {

        self.computed_polygon = Some(self.create_polygon());
        self
       
    }

    pub fn parse_geometry(&self, coord_string: &String) -> Vec<Coord<f64>> {
        let coordinates_str: Vec<&str> = coord_string.split(" ").collect();

        // Parse coordinates into a Vec of `Coord<f64>`
        let mut coords: Vec<Coord<f64>> = Vec::new();

        for coordinate in coordinates_str {
            let parts: Vec<&str> = coordinate.split(',').collect();
            if parts.len() == 2 {
                let e: f64 = parts[0].parse().expect("Invalid x coordinate");
                let n: f64 = parts[1].parse().expect("Invalid y coordinate");

                let (lon, lat) = self.proj.transform(e, n);
                coords.push(Coord { x: lon, y: lat });
            } else {
                println!("Invalid coordinate format: {}", coordinate);
            }
        }

        coords
    }
 
    pub fn get_geometries(&self) -> (LineString, Vec<LineString>) {
        let polygon = &self
            .stand_basic_data
            .polygon_geometry
            .polygon_property
            .polygon;

        let exterior = &polygon.exterior.linear_ring.coordinates;
        let exterior_geometry = LineString::new(self.parse_geometry(&exterior).to_owned());

        let interior_geometry: Vec<LineString> = polygon
            .interior
            .iter()
            .map(|f| {
                let geometry = self.parse_geometry(&f.linear_ring.coordinates);
                LineString::new(geometry)
            })
            .collect();

        (exterior_geometry, interior_geometry)
    }

    pub fn create_polygon(&mut self) -> Polygon {
        // Projection from ETRS-TM35FIN to WGS84
        self.proj = Projection::new(CRS::Epsg3067, CRS::Epsg4326);

        let (exterior, interior) = self.get_geometries();

        let polygon = Polygon::new(exterior, interior);

        polygon
    }

    pub fn summary_stem_count(&self) -> Option<u32> {

        let last_data_date = match self.get_last_tree_stand_data_date() {
            Some(data) => data,
            None => return None 
        };

        match &last_data_date.tree_stand_summary {
            Some(summary) => return Some(summary.stem_count),
            None => return None
        };
        
    }

    pub fn stem_count_in_stratum(&self) -> bool {
        let stratums = self.get_stratums();

        let stratum_vec = match stratums {
            Some(stratum) => stratum,
            None => return false
        };

        true
    }

    pub fn get_stratums(&self) -> Option<Vec<TreeStratum>> {
        let last_data_date = match self.get_last_tree_stand_data_date() {
            Some(data) => data,
            None => return None 
        };

        let stratums = last_data_date.tree_strata.tree_stratum.to_owned();
        Some(stratums)
    }

    pub fn get_strata(&self) -> Option<TreeStrata> {
        let last_data_date = match self.get_last_tree_stand_data_date() {
            Some(data) => data,
            None => return None 
        };

        let strata = &last_data_date.tree_strata.tree_stratum;
        let strata = TreeStrata::new(strata.to_vec());
        Some(strata)
    }

    pub fn get_last_tree_stand_data_date(&self) -> Option<TreeStandDataDate> {
        let stand_data = match &self.tree_stand_data {
            Some(data) => data,
            None => return None 
        };

        match stand_data.tree_stand_data_date.last() {
            Some(last_data_date) => Some(last_data_date.to_owned()),
            None => return None 
        }
    }
}
