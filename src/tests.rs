#[test]
fn test_writing_to_json() {
    let test_json_path = "test_json_from_xml.json";

    let xml_property = ForestPropertyData::from_xml_file("forestpropertydata.xml");

    match xml_property.write_to_json_file(&test_json_path) {
        std::result::Result::Ok(_) => assert!(true),
        _ => assert!(false),
    }

    fs::remove_file(test_json_path).unwrap()
}

#[test]
fn test_parsers() {
    let xml_property = ForestPropertyData::from_xml_file("forestpropertydata.xml");
    xml_property
        .write_to_json_file("json_from_xml.json")
        .expect("writing JSON failed");

    let json_property = ForestPropertyData::from_json_file("json_from_xml.json");

    let xml_real_estate = xml_property.real_estates.real_estate.first().unwrap();
    let json_real_estate = json_property.real_estates.real_estate.first().unwrap();

    let xml_id = &xml_real_estate.id;
    let json_id = &json_real_estate.id;

    assert!(xml_id == json_id, "JSON and XML file parsing");

    let xml_stands = xml_real_estate.get_stands();
    let json_stands = json_real_estate.get_stands();

    for i in 0..xml_stands.iter().len() {
        assert!(
            xml_stands[i].id == json_stands[i].id,
            "stand is matches with id: {}",
            i
        )
    }
}

// Run wih `cargo test -- --nocapture` to see the print statements
#[test]
fn test_find_stands_in_bounding_box() {
    let property = ForestPropertyData::from_xml_file("forestpropertydata.xml");
    let real_estate = property.real_estates.real_estate[0].clone();
    let all_stands = real_estate.get_stands();

    let mut stands = Vec::new();
    for stand in all_stands {
        stands.push(stand.clone());
    }
    println!("\nTotal stands: {:?}", stands.len());
    let min_x = 425400.0;
    let max_x = min_x + 6000.0;
    let min_y = 7369000.0;
    let max_y = min_y + 6000.0;

    let bbox = geo::Polygon::new(
        LineString(vec![
            Coord { x: min_x, y: min_y },
            Coord { x: max_x, y: min_y },
            Coord { x: max_x, y: max_y },
            Coord { x: min_x, y: max_y },
            Coord { x: min_x, y: min_y },
        ]),
        vec![],
    );
    let stands = find_stands_in_bounding_box(&stands, &bbox);

    if !stands.is_none() {
        println!(
            "Stands in bounding box: {:?}",
            stands.clone().unwrap().len()
        );
        for stand in &stands.unwrap() {
            println!("Stand number {:?}", stand.stand_basic_data.stand_number);
        }
    }
}
