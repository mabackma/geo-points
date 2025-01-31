use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct PolygonGeometry {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "pointProperty")]
    pub point_property: PointProperty,
    #[serde(rename = "polygonProperty")]
    pub polygon_property: PolygonProperty,
}

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct PointProperty {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "Point")]
    pub point: Point,
}

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Point {
    #[serde(rename = "@srsName", default)]
    pub srs_name: String,
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "coordinates")]
    pub coordinates: String,
}

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct PolygonProperty {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "Polygon")]
    pub polygon: Polygon,
}

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Polygon {
    #[serde(rename = "@srsName", default)]
    pub srs_name: String,
    #[serde(rename = "$text", default)]
    pub text: Option<String>,
    #[serde(rename = "interior", default)]
    pub interior: Vec<Interior>,
    #[serde(rename = "exterior")]
    pub exterior: Exterior,
}

impl Polygon {
    pub fn to_geo_polygon(&self) -> geo_types::Polygon<f64> {
        let exterior = self.exterior.linear_ring.to_geo_line_string();
        let interiors = self.interior.iter().map(|i| i.linear_ring.to_geo_line_string()).collect();
        
        geo_types::Polygon::new(exterior, interiors)
    }
}

#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Interior {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "LinearRing")]
    pub linear_ring: LinearRing,
}


#[derive(Serialize, Deserialize, Debug,Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct Exterior {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "LinearRing")]
    pub linear_ring: LinearRing,
}

#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct LinearRing {
    #[serde(rename = "$text")]
    pub text: Option<String>,
    #[serde(rename = "coordinates")]
    pub coordinates: String,
}

use geo_types::{LineString, Coordinate};

impl LinearRing {
    pub fn to_geo_line_string(&self) -> LineString<f64> {
        let coords: Vec<Coordinate<f64>> = self.coordinates.split_whitespace().map(|c| {
            let mut split = c.split(",");
            let x = split.next().unwrap().parse::<f64>().unwrap();
            let y = split.next().unwrap().parse::<f64>().unwrap();
            Coordinate { x, y } // Create a Coord from x and y
        }).collect();

        LineString(coords) // Create a LineString from the vector of Coord
    }
}
