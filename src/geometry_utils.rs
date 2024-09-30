use crate::forest_property::tree_stand_data::TreeStrata;
use crate::forest_property::tree::Tree;
use crate::jittered_hexagonal_sampling::{GridOptions, JitteredHexagonalGridSampling};
use crate::projection::{Projection, CRS};

use geo_types::Polygon;
use geo::{BoundingRect, Coord, LineString};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use core::f32::consts::PI;

// Get minimum and maximum x and y coordinates of a polygon
pub fn get_min_max_coordinates(p: &Polygon<f64>) -> (f64, f64, f64, f64) {
    let rect = p.bounding_rect().unwrap();
    let min_x = rect.min().x;
    let max_x = rect.max().x;
    let min_y = rect.min().y;
    let max_y = rect.max().y;

    (min_x, max_x, min_y, max_y)
}

pub fn generate_radius(total_stem_count: u32, area: f32) -> f32 {
    let total_trees = total_stem_count as f32 * area / 10000.0;

    let mut ratio_fix = 1.3;

    if total_trees < 250.0 {
        ratio_fix = ((total_trees / 250.0) * 0.6) + 1.3;
    }
    let square_to_circle_ratio = 1.273 / ratio_fix;

    let tree_needed_area = area / total_trees / square_to_circle_ratio;
    // Calculate the radius based on the mean height of the tree species
    (tree_needed_area / PI).sqrt()
}

// Generates random trees for all strata with jittered grid sampling
pub fn generate_random_trees(p: &Polygon, strata: &TreeStrata, area_ratio: f64) -> Vec<Tree> {
    let total_stem_count = strata.tree_stratum.iter().fold(0, |mut acc: u32, f| {
        acc += f.stem_count;
        acc
    });

    let trees = strata
        .tree_stratum
        .par_iter()
        .map(|stratum| {
            let tree_amount = (stratum.stem_count as f64) * area_ratio;
            let amount = tree_amount.round() as u32;

            let mut radius = generate_radius(
                total_stem_count,
                stratum.basal_area
            );
            
            radius *= 0.00001;

            // Jittered Grid Version 2
            let rng = rand::thread_rng();
            let options = GridOptions {
                polygon: p.to_owned(),
                radius: (radius).into(),
                jitter: Some(0.6666),
                point_limit: Some(amount as usize),
            };

            let mut grid = JitteredHexagonalGridSampling::new(rng, options);

            let points =  grid.fill();

            if points.len() == 0 {
                //println!("\tNo trees generated for stratum with basal area {}, stem count {}, mean height {}", stratum.basal_area, stratum.stem_count, stratum.mean_height);
            }
            else if points.len() < amount as usize {
                println!("Generated {} / {} trees for stratum with basal area {}, stem count {}, mean height {}.", points.len(), amount, stratum.basal_area, stratum.stem_count, stratum.mean_height);
            }

            let trees_strata: Vec<Tree> = points.iter().map(|pair: &[f64; 2]| {
                Tree::new(
                    stratum.tree_species,
                    stratum.mean_height,
                    (pair[0], pair[1], 0.0),
                )
            }).collect();
            trees_strata
        })
        .flatten();

    trees.collect()
}

pub fn polygon_to_wgs84(p: &Polygon) -> Polygon {
    let proj = Projection::new(CRS::Epsg3067, CRS::Epsg4326);
    let mut coords: Vec<Coord<f64>> = Vec::new();

    for coordinate in p.exterior().points() {
        let e: f64 = coordinate.x();
        let n: f64 = coordinate.y();
        
        let (lon, lat) = proj.transform(e, n);
        coords.push(Coord { x: lon, y: lat });
    }

    Polygon::new(LineString::from(coords), vec![])
}