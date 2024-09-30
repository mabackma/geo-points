use geo::{algorithm::contains::Contains, prelude::*, Polygon, Coord};
use rand::seq::SliceRandom;
use rand::Rng;

const HEX_SIDE: f64 = 0.8660254037844386;

const VERTICES: [[f64; 2]; 6] = [
    [0.0, -1.0],
    [HEX_SIDE, -0.5],
    [HEX_SIDE, 0.5],
    [0.0, 1.0],
    [-HEX_SIDE, 0.5],
    [-HEX_SIDE, -0.5],
];

pub struct JitteredHexagonalGridSampling<R: Rng> {
    polygon: Polygon<f64>,
    r: f64,
    jitter: f64,
    jitter_radius: f64,
    rng: R,
    max_y: usize,
    max_x_even: usize,
    max_x_odd: usize,
    current_x: usize,
    current_y: usize,
    min_x: f64,
    min_y: f64,
    sample_points: Vec<[f64; 2]>,
    point_limit: Option<usize>,
}

impl<R: Rng> JitteredHexagonalGridSampling<R> {
    pub fn new(rng: R, options: GridOptions) -> Self {
        let r = options.radius;
        let jitter = options.jitter.unwrap_or(0.6666);
        let jitter_radius = r * jitter;

        let bounding_rect = options.polygon.bounding_rect().unwrap();
        let min_x = bounding_rect.min().x;
        let min_y = bounding_rect.min().y;

        let max_y = ((bounding_rect.max().y - min_y) / r).ceil() as usize;
        let max_x_even = ((bounding_rect.max().x - min_x) / (r * 2.0 * HEX_SIDE) + 0.5).ceil() as usize;
        let max_x_odd = ((bounding_rect.max().x - min_x) / (r * 2.0 * HEX_SIDE)).ceil() as usize;

        Self {
            polygon: options.polygon,
            r,
            jitter,
            jitter_radius,
            rng,
            max_y,
            max_x_even,
            max_x_odd,
            current_x: 0,
            current_y: 0,
            min_x,
            min_y,
            sample_points: Vec::new(),
            point_limit: options.point_limit,
        }
    }

    pub fn get_all_points(&self) -> &[[f64; 2]] {
        &self.sample_points
    }
 
    pub fn generate_all_points(&mut self) {
        while self.current_y < self.max_y {
            let y_odd = self.current_y % 2 == 1;
            let max_x = if y_odd { self.max_x_odd } else { self.max_x_even };
            while self.current_x < max_x {
                let cx = self.min_x + (self.current_x as f64 * 2.0 + if y_odd { 1.0 } else { 0.0 }) * self.r * HEX_SIDE;
                let cy = self.min_y + self.current_y as f64 * 1.5 * self.r;

                let mut p = self.rng.gen_range(0.0..6.0);
                let q = self.rng.gen::<f64>();
                let tri = p as usize;
                p %= 1.0;

                let v0 = VERTICES[tri];
                let v1 = VERTICES[(tri + 1) % 6];

                let (p, q) = if p + q > 1.0 { (1.0 - p, 1.0 - q) } else { (p, q) };

                let point = [
                    cx + (v0[0] * p + v1[0] * q) * self.jitter_radius,
                    cy + (v0[1] * p + v1[1] * q) * self.jitter_radius,
                ];

                if self.polygon.contains(&Coord { x: point[0], y: point[1] }) {
                    self.sample_points.push(point);
                }

                self.current_x += 1;
            }

            self.current_x = 0;
            self.current_y += 1;
        }
    }

    pub fn fill(&mut self) -> Vec<[f64; 2]> {
 
        self.generate_all_points();


        if let Some(limit) = self.point_limit {
            let mut sampled_points = self.sample_points.clone();
            sampled_points.shuffle(&mut self.rng);
            sampled_points.truncate(limit);
            sampled_points
        } else {
            self.sample_points.clone()
        }
    }

}

pub struct GridOptions {
    pub polygon: Polygon<f64>,
    pub radius: f64,
    pub jitter: Option<f64>,
    pub point_limit: Option<usize>,
}