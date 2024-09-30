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