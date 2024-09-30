use std::fmt;
use proj4rs::proj::Proj;

pub struct Projection {
    from: Proj,
    to: Proj
}

pub const EPSG_3067: &'static str  = "+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs";
pub const EPSG_4326: &'static str  = "+proj=longlat +datum=WGS84 +no_defs +type=crs";

pub enum CRS {
    Epsg3067,
    Epsg4326
}

impl CRS {
    
    fn proj_str(&self) -> &str {

        match &self {
            CRS::Epsg3067 => &EPSG_3067,
            CRS::Epsg4326 => &EPSG_4326
        }

    }
}

impl Projection {

    pub fn new(from: CRS, to: CRS) -> Self {

        Projection {
            from: Proj::from_proj_string(from.proj_str()).unwrap(),
            to: Proj::from_proj_string(to.proj_str()).unwrap()
        }
    }

    pub fn set_projection_from(&mut self, crs: CRS) -> &Self {

        self.from = Proj::from_proj_string(crs.proj_str()).unwrap();
        self
    }

    pub fn set_projection_to(&mut self, crs: CRS) -> &Self {
        self.to = Proj::from_proj_string(crs.proj_str()).unwrap();
        self
    }

    pub fn transform(&self, x:f64, y:f64) -> (f64, f64) {

        let mut point_3d = (x, y, 0.0);

        proj4rs::transform::transform(&self.from, &self.to, &mut point_3d).unwrap();
    
        // Note that angular unit is radians, not degrees
        (point_3d.0.to_degrees(), point_3d.1.to_degrees())
    }

    pub fn transform_back(&self, x:f64, y:f64) -> (f64, f64) {
        // Convert degrees to radians
        let mut point_3d = (x.to_radians(), y.to_radians(), 0.0);

        proj4rs::transform::transform(&self.from, &self.to, &mut point_3d).unwrap();

        (point_3d.0, point_3d.1)
    }
}

// Implement Debug for Projection
impl fmt::Debug for Projection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Projection {{ from: {:?}, to: {:?} }}", self.from, self.to)
    }
}

impl PartialEq for Projection {
    fn eq(&self, other: &Self) -> bool {
        //self.from == other.from && self.to == other.to
        true
    }
}

// Implement Clone for Projection
impl Clone for Projection {
    fn clone(&self) -> Self {
        Projection {
            from: self.from.clone(),
            to: self.to.clone(),
        }
    }
}

impl Default for Projection {
    fn default() -> Self {
        Projection {
            from: Proj::from_proj_string("+proj=latlong").unwrap(),
            to: Proj::from_proj_string("+proj=latlong").unwrap(),
        }
    }
}

#[test]
fn test_projection() {
    // EPSG:3067 - TM35FIN(E,N) -- Finland
    let from = Proj::from_proj_string(
        "+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs",
    )
    .unwrap();

    // EPSG:4326 - WGS84
    let to = Proj::from_proj_string("+proj=longlat +datum=WGS84 +no_defs +type=crs")
    .unwrap();
    
    /*  N=7369564.333, E=427997.035 */
    let epsg3067_northern = 7369564.333;
    let epsg3067_eastern = 427997.035;


    let mut point_3d = (epsg3067_eastern, epsg3067_northern, 0.0);
    proj4rs::transform::transform(&from, &to, &mut point_3d).unwrap();

    // Note that angular unit is radians, not degrees
    let (longitude, latitude,_height) = point_3d;

    // Output in longitude, latitude
    println!("LatLng: {},{}", latitude.to_degrees(), longitude.to_degrees());

    // Projection validated here:
    // https://epsg.io/transform#s_srs=3067&t_srs=4326&ops=1149&x=427997.0350000&y=7369564.3330000
    assert!(((latitude.to_degrees() * 1E6).round() / 1E6) == 66.437124, "projection match with lat");
    assert!(((longitude.to_degrees() * 1E6).round() / 1E6) == 25.385742, "projection match with lon");

    // EPSG:3067 E 427997.035 -> EPSG:4326 longitude 25°23'8.67" = 25.385742
    // EPSG:3067 N 7369564.333 -> EPSG:4326 latitude 66°26'13.646" = 66.437124
}

#[test]
fn test_projection_impl() {
    // EPSG:3067 - TM35FIN(E,N) -- Finland
    let proj = Projection::new(CRS::Epsg3067, CRS::Epsg3067);

    
    /*  N=7369564.333, E=427997.035 */
    let epsg3067_northern = 7369564.333;
    let epsg3067_eastern = 427997.035;

    let (lon, lat) = proj.transform(epsg3067_northern, epsg3067_eastern);

    // Output in longitude, latitude
    println!("LatLng: {},{}", lat, lon);

    // Projection validated here:
    // https://epsg.io/transform#s_srs=3067&t_srs=4326&ops=1149&x=427997.0350000&y=7369564.3330000
    assert!(((lat * 1E6).round() / 1E6) == 66.437124, "projection match with lat");
    assert!(((lon * 1E6).round() / 1E6) == 25.385742, "projection match with lon");

    // EPSG:3067 E 427997.035 -> EPSG:4326 longitude 25°23'8.67" = 25.385742
    // EPSG:3067 N 7369564.333 -> EPSG:4326 latitude 66°26'13.646" = 66.437124
}