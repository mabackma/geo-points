use serde_json::Value;
use serde::{Deserialize, Serialize};

use crate::forest_property::forest_property_data::TreeStratum;

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TreeStandData {
    #[serde(rename = "TreeStandDataDate")]
    pub tree_stand_data_date: Vec<TreeStandDataDate>,
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TreeStandDataDate {
    #[serde(rename = "TreeStrata")]
    pub tree_strata: TreeStrata,
    #[serde(rename = "TreeStandSummary")]
    pub tree_stand_summary: Option<TreeStandSummary>,
    #[serde(rename = "DeadTreeStrata")]
    pub dead_tree_strata: Option<DeadTreeStrata>,
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TreeStrata {
    #[serde(rename = "TreeStratum")]
    pub tree_stratum: Vec<TreeStratum>,
}

impl TreeStrata {
    pub fn new(tree_stratum: Vec<TreeStratum>) -> Self {
        TreeStrata { tree_stratum }
    }
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TreeStandSummary {
    #[serde(rename = "ChangeState")]
    pub change_state: i64,
    #[serde(rename = "MeanAge")]
    pub mean_age: i64,
    #[serde(rename = "BasalArea")]
    pub basal_area: f64,
    #[serde(rename = "StemCount")]
    pub stem_count: i64,
    #[serde(rename = "MeanDiameter")]
    pub mean_diameter: f64,
    #[serde(rename = "MeanHeight")]
    pub mean_height: f64,
    #[serde(rename = "Volume")]
    pub volume: f64,
    #[serde(rename = "VolumeGrowth")]
    pub volume_growth: f64,
    #[serde(rename = "Value")]
    pub value: Option<f64>,
    #[serde(rename = "ValueGrowthPercent")]
    pub value_growth_percent: Option<f64>,
    #[serde(rename = "SawLogVolume")]
    pub saw_log_volume: Option<f64>,
    #[serde(rename = "PulpWoodVolume")]
    pub pulp_wood_volume: Option<f64>,
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DeadTreeStrata {
    #[serde(rename = "DeadTreeStratum")]
    pub dead_tree_stratum: Value,
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Operations {
    #[serde(rename = "Operation")]
    pub operation: Value,
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpecialFeatures {
    #[serde(rename = "SpecialFeature")]
    pub special_feature: Value,
}