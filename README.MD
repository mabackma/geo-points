## Building and Running the Project

1. To build the project:

    ```bash
    cargo build --release
    ```

2. To run the executable:

    - On Windows:
        ```bash
        target/release/geo_points.exe
        ```

    - On macOS/Linux:
        ```bash
        ./target/release/geo_points
        ```

---

The executable will generate a random bounding box and draw the stands and the generated trees inside the bounding box into stands_in_bbox_image.png. The user will then be asked to choose a stand from the XML file. Choose 0 for Realestate and then choose the stand number to draw the stand into selected_stand_image.png.