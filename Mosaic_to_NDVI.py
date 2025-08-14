import os
import re
import rasterio
import numpy
from xml.dom import minidom

# Input folder containing GeoTIFF files (Planet outputs)
input_folder = '/Users/u0150402/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Output_planet_mosaic'

# Output folder for NDVI files
output_folder = '/Volumes/T7 Shield/NDVI_mosaic_output'

# Find all GeoTIFF files starting with "new_mosaic_by_date_" and ending with ".tif"
tif_files = [f for f in os.listdir(input_folder) if f.startswith("new_mosaic_by_date_") and f.endswith(".tif")]

for tif_file in tif_files:
    # Construct the full file path for the input image
    image_file = os.path.join(input_folder, tif_file)
    
    # Extract date and code from the GeoTIFF file name (adjust as per your naming convention)
    date_code = os.path.splitext(tif_file)[0]

    # Load red and NIR bands (assuming bands are in positions 3 and 4)
    with rasterio.open(image_file) as src:
        band_red = src.read(3)
        band_nir = src.read(4)

    # Calculate NDVI
    ndvi = (band_nir.astype(float) - band_red.astype(float)) / (band_nir + band_red)

    # Set spatial characteristics of the output to mirror the input
    kwargs = src.meta
    kwargs.update(
        dtype=rasterio.float32,
        count=1
    )

    # Specify the output path using the date and code
    output_path = os.path.join(output_folder, f"{date_code}_NDVI.tif")

    # Write the NDVI file
    with rasterio.open(output_path, 'w', **kwargs) as dst:
        dst.write_band(1, ndvi.astype(rasterio.float32))

    print(f"NDVI file for {date_code} saved at: {output_path}")
