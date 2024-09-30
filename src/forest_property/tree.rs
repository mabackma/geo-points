#[derive(Default, Debug, Clone, Copy)]
pub struct Tree {
    species: u8,
    mean_height: f32,
    position: (f64, f64, f64),
}

impl Tree {
    pub fn new(species: u8, mean_height: f32, position: (f64, f64, f64)) -> Self {
        Tree {
            species,
            mean_height,
            position,
        }
    }

    pub fn species(&self) -> u8 {
        self.species
    }

    pub fn position(&self) -> (f64, f64, f64) {
        self.position
    }
}