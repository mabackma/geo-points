## Building and Running the Project

1. To build the project:

    ```bash
    cargo build --release
    ```

2. To run the executable:

    ```bash
    ./target/release/geo_points
    ```
---

The executable will generate a random bounding box and save the stands and the generated trees inside the bounding box into an image file named stands_in_bbox_image.png. It will also save them as geojson as stands_in_bbox.geojson.

The user will then be asked to choose a stand from the XML file. Choose 0 for Realestate and then choose the stand number. It will save the stand into an image file named selected_stand_image.png and its geojson as selected_stand.geojson.