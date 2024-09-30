use std::fs::{self, File};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use super::{geometry::PolygonGeometry, stand::{Stand, Stands}};

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ForestPropertyData {
    #[serde(rename = "RealEstates")]
    pub real_estates: RealEstates,
}

pub fn read_number_cli(min: usize, max: usize) -> usize {
    let mut buf = String::new();

    std::io::stdin()
        .read_line(&mut buf)
        .expect("Failed to read line");

    let number: usize = buf.trim().parse().expect("Please type a number!");

    if number > max || number < min {
        panic!("Chosen number out of range")
    }

    number
}

impl ForestPropertyData {
    pub fn from_xml_file(path: &str) -> ForestPropertyData {
        let xml = fs::read_to_string(path).expect("Could not read the XML file");
        ForestPropertyData::parse_from_str(xml.as_str())
    }

    pub fn from_xml_str(xml_str: &str) -> ForestPropertyData {
        ForestPropertyData::parse_from_str(xml_str)
    }

    #[cfg(test)]
    pub fn from_json_file(path: &str) -> ForestPropertyData {
        let json = fs::read_to_string(path).expect("Could not read the JSON file");

        
        let property : ForestPropertyData = serde_json::from_str(json.as_str()).expect("Error parsing json file");
        
        property
    }

    #[cfg(test)]
    pub fn write_to_json_file(&self, path: &str) -> anyhow::Result<(), anyhow::Error> {
        let json_file = File::create(path)?;

        serde_json::to_writer(json_file, &self)?;

        Ok(())
    }

    pub fn parse_from_str(xml: &str) -> ForestPropertyData {
        quick_xml::de::from_str(xml).expect("Could not parse the XML")
    }

    // Parcels are not probably needed in this context but its good to keep them just in case
    pub fn choose_parcel(&self) -> Parcel {
        let parcels: &Vec<Parcel> = &self.real_estates.real_estate.first().unwrap().parcels.parcel;

        println!("\nParcels:");
        for (i, parcel) in parcels.iter().enumerate() {
            print!("{}. {:?}, ", i, parcel.parcel_number);
        }
        println!("Choose a parcel to view: ");

        let parcel_index = read_number_cli(0, parcels.len());
        parcels[parcel_index].to_owned()
    }

    pub fn get_stand_cli(&self) -> Stand {
        let real_estates = &self.real_estates.real_estate;

        println!("Realestates:");
        for (i, real_estate) in real_estates.iter().enumerate() {
            println!("{}. {:?}, ", i.to_string(), real_estate.real_estate_name);
        }
        println!("Choose a realestate number to view: ");

        let real_estate_index = read_number_cli(0, real_estates.len());

        let real_estate = &self.real_estates.real_estate[real_estate_index];

        let stands_data: Vec<Stand> = real_estate.get_stands();

        for (i, stand) in stands_data.iter().enumerate() {
            let StandBasicData {
                area, stand_number, ..
            } = &stand.stand_basic_data;

            println!(
                "{}. {:?}, {}ha, stand number: {}",
                i, stand.id, area, stand_number
            )
        }

        println!("Choose a stand:");

        let stand_index = read_number_cli(0, stands_data.len());

        let stand = &stands_data[stand_index];

        stand.to_owned()
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RealEstates {
    
    #[serde(rename = "RealEstate")]
    pub real_estate: Vec<RealEstate>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RealEstate {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "MunicipalityNumber")]
    pub municipality_number: u16,
    #[serde(rename = "AreaNumber")]
    pub area_number: u16,
    #[serde(rename = "GroupNumber")]
    pub group_number: u16,
    #[serde(rename = "UnitNumber")]
    pub unit_number: u16,
    #[serde(rename = "RealEstateName")]
    pub real_estate_name: String,
    #[serde(rename = "Parcels")]
    pub parcels: Parcels,
}

impl RealEstate {
    
    pub fn get_stands(&self) -> Vec<Stand> {

        let parcels = &self.parcels.parcel;

        let stands_data: Vec<Stand> = parcels
            .into_par_iter()
            .map(|parcel: &Parcel| {
                let stands: Vec<Stand> = parcel.stands.stand.iter().map( | f| f.to_owned().compute_polygon().to_owned()).collect();
                stands
            }).flatten()
            .collect();

        stands_data
    } 

}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Parcels {
    
    #[serde(rename = "Parcel")]
    pub parcel: Vec<Parcel>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Parcel {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "ParcelNumber")]
    pub parcel_number: i64,
    #[serde(rename = "Stands")]
    pub stands: Stands,
}


#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct StandBasicData {
    
    #[serde(rename = "Identifiers")]
    pub identifiers: Option<Identifiers>,
    #[serde(rename = "CuttingRestriction", default)]
    pub cutting_restriction: u8,
    #[serde(rename = "StandInfo")]
    pub stand_info: Option<String>,
    #[serde(rename = "DitchingYear")]
    pub ditching_year: Option<u16>,
    #[serde(rename = "ChangeTime")]
    pub change_time: String,
    #[serde(rename = "CompleteState")]
    pub complete_state: u8,
    #[serde(rename = "StandNumber")]
    pub stand_number: u16,
    #[serde(rename = "StandNumberExtension")]
    pub stand_number_extension: String,
    #[serde(rename = "MainGroup")]
    pub main_group: u8,
    #[serde(rename = "StandBasicDataDate")]
    pub stand_basic_data_date: String,
    #[serde(rename = "Area")]
    pub area: f32,
    #[serde(rename = "PolygonGeometry")]
    pub polygon_geometry: PolygonGeometry,
    #[serde(rename = "AreaDecrease")]
    pub area_decrease: Option<String>,
    #[serde(rename = "Accessibility")]
    pub accessibility: Option<u8>,
    #[serde(rename = "MainTreeSpecies")]
    pub main_tree_species: Option<u8>,
    #[serde(rename = "StandQuality")]
    pub stand_quality: Option<u8>,
    #[serde(rename = "DevelopmentClass")]
    pub development_class: Option<String>,
    #[serde(rename = "DrainageState")]
    pub drainage_state: Option<u8>,
    #[serde(rename = "SoilType")]
    pub soil_type: Option<u8>,
    #[serde(rename = "FertilityClass")]
    pub fertility_class: Option<u8>,
    #[serde(rename = "SubGroup")]
    pub sub_group: Option<u8>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Identifiers {
    
    #[serde(rename = "Identifier")]
    pub identifier: Identifier,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Identifier {
    
    #[serde(rename = "IdentifierType")]
    pub identifier_type: String,
    #[serde(rename = "IdentifierValue")]
    pub identifier_value: String,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SpecialFeatures {
    
    #[serde(rename = "SpecialFeature")]
    pub special_feature: Vec<SpecialFeature>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SpecialFeature {
    #[serde(rename = "@id")]
    pub id: u32,
    
    #[serde(rename = "FeatureAdditionalCode")]
    pub feature_additional_code: Option<String>,
    #[serde(rename = "FeatureCode")]
    pub feature_code: String,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Operations {
    
    #[serde(rename = "Operation")]
    pub operation: Vec<Operation>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Operation {
    #[serde(rename = "@mainType")]
    pub main_type: String,
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "CompletionData")]
    pub completion_data: Option<CompletionData>,
    #[serde(rename = "OperationInfo")]
    pub operation_info: Option<String>,
    #[serde(rename = "Specifications")]
    pub specifications: Option<Specifications>,
    #[serde(rename = "Silviculture")]
    pub silviculture: Option<Silviculture>,
    #[serde(rename = "ChangeTime")]
    pub change_time: String,
    #[serde(rename = "OperationType")]
    pub operation_type: String,
    #[serde(rename = "ProposalData")]
    pub proposal_data: ProposalData,
    #[serde(rename = "Cutting")]
    pub cutting: Option<Cutting>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct CompletionData {
    
    #[serde(rename = "CompletionDate")]
    pub completion_date: String,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Specifications {
    
    #[serde(rename = "Specification")]
    pub specification: Vec<Specification>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Specification {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "SpecificationCode")]
    pub specification_code: u32,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Silviculture {}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ProposalData {
    #[serde(rename = "ProposalType")]
    pub proposal_type: u32,
    #[serde(rename = "ProposalYear")]
    pub proposal_year: u32,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Cutting {
    
    #[serde(rename = "CuttingVolume",default = "default_zero_f32")]
    pub cutting_volume: f32,
    #[serde(rename = "Assortments")]
    pub assortments: Option<Assortments>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Assortments {
    
    #[serde(rename = "Assortment", default)]
    pub assortment: Vec<Assortment>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Assortment {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "TreeSpecies")]
    pub tree_species: Option<String>, // This fails to parse even with default values unless type is Option<String>. Requires manual parsing.
    #[serde(rename = "StemType", default = "default_zero_u32")]
    pub stem_type: u32,
    #[serde(rename = "AssortmentVolume", default = "default_zero_f32")]
    pub assortment_volume: f32,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct TreeStandData {
    
    #[serde(rename = "TreeStandDataDate")]
    pub tree_stand_data_date: Vec<TreeStandDataDate>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct TreeStandDataDate {
    #[serde(rename = "@date", default)]
    pub date: String,
    #[serde(rename = "@type", default)]
    pub tree_stand_data_date_type: u8,
    
    #[serde(rename = "DeadTreeStrata")]
    pub dead_tree_strata: Option<DeadTreeStrata>,
    #[serde(rename = "TreeStrata")]
    pub tree_strata: TreeStrata,
    #[serde(rename = "TreeStandSummary")]
    pub tree_stand_summary: Option<TreeStandSummary>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct DeadTreeStrata {
    
    #[serde(rename = "DeadTreeStratum")]
    pub dead_tree_stratum: Vec<DeadTreeStratum>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Copy)]
pub struct DeadTreeStratum {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "DeadTreeType")]
    pub dead_tree_type: u8,
    #[serde(rename = "TreeSpecies")]
    pub tree_species: u8,
    #[serde(rename = "Volume", default = "default_zero_f32")]
    pub volume: f32,
    #[serde(rename = "MeanDiameter", default = "default_zero_f32")]
    pub mean_diameter: f32,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct TreeStrata {
    
    #[serde(rename = "TreeStratum")]
    pub tree_stratum: Vec<TreeStratum>,
}


fn default_zero_u32() -> u32 {
    0
}

fn default_zero_f32() -> f32 {
    0.0
}

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize, Copy)]
#[serde(rename_all = "camelCase")]
pub struct TreeStratum {
    #[serde(rename = "@id", default)]
    pub id: u32,
    #[serde(rename = "StratumNumber")]
    pub stratum_number: u32,
    #[serde(rename = "TreeSpecies")]
    pub tree_species: u8,
    #[serde(rename = "Storey")]
    pub storey: u8,
    #[serde(rename = "Age")]
    pub age: u8,
    #[serde(rename = "StemCount", default = "default_zero_u32")]
    pub stem_count: u32,
    #[serde(rename = "MeanDiameter", default = "default_zero_f32")]
    pub mean_diameter: f32,
    #[serde(rename = "MeanHeight")]
    pub mean_height: f32,
    #[serde(rename = "DataSource")]
    pub data_source: u32,
    #[serde(rename = "BasalArea",  default = "default_zero_f32")]
    pub basal_area: f32,
    #[serde(rename = "SawLogPercent", default = "default_zero_f32")]
    pub saw_log_percent: f32,
    #[serde(rename = "SawLogVolume", default = "default_zero_f32")]
    pub saw_log_volume: f32,
    #[serde(rename = "VolumeGrowth", default = "default_zero_f32")]
    pub volume_growth: f32,
    #[serde(rename = "Volume", default = "default_zero_f32")]
    pub volume: f32,
    #[serde(rename = "PulpWoodVolume", default = "default_zero_f32")]
    pub pulp_wood_volume: f32,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Copy)]
pub struct TreeStandSummary {
    #[serde(rename = "@id")]
    pub id: u32,
    #[serde(rename = "PulpWoodVolume", default = "default_zero_f32")]
    pub pulp_wood_volume: f32,
    #[serde(rename = "SawLogVolume", default = "default_zero_f32")]
    pub saw_log_volume: f32,
    #[serde(rename = "MeanAge")]
    pub mean_age: f32,
    #[serde(rename = "BasalArea")]
    pub basal_area: f32,
    #[serde(rename = "StemCount")]
    pub stem_count: u32,
    #[serde(rename = "MeanDiameter")]
    pub mean_diameter: f32,
    #[serde(rename = "MeanHeight")]
    pub mean_height: f32,
    #[serde(rename = "Volume")]
    pub volume: f32,
    #[serde(rename = "VolumeGrowth")]
    pub volume_growth: f32,
    #[serde(rename = "ValueGrowthPercent", default = "default_zero_f32")]
    pub value_growth_percent: f32,
}