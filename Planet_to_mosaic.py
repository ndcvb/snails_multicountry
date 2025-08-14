import os
import re
import numpy as np
from os.path import join
from rasterio.merge import merge
from rasterio import open as rasterio_open

# Input and output folders
input_folder = '/Users/u0150402/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Output_planet'
output_folder = '/Users/u0150402/Library/CloudStorage/OneDrive-KULeuven/Uganda_Congo/Data/Uganda/Output_planet_mosaic'

# Create the output folder if it does not exist
os.makedirs(output_folder, exist_ok=True)

# Find all GeoTIFF files with the ending "AnalyticMS_SR_harmonized_clip.tif"
tif_files = [f for f in os.listdir(input_folder) if f.endswith('AnalyticMS_SR_harmonized_clip.tif')]

# Extract dates from file names and group files by date
date_file_dict = {}

# Updated regex pattern
date_pattern = re.compile(r"(\d{8})_\d{6}")

for tif_file in tif_files:
    # Construct the full file paths
    image_file = join(input_folder, tif_file)
    
    # Extract the date part from the filename
    match = date_pattern.search(tif_file)
    if match:
        date_str = match.group(1)
        
        # Group by date
        if date_str not in date_file_dict:
            date_file_dict[date_str] = []
        date_file_dict[date_str].append(image_file)
    else:
        print(f"No date found in file: {tif_file}")

# Process each date group
for idx, (date, files) in enumerate(date_file_dict.items(), start=1):
    print(f"Processing Mosaic {idx} for date: {date}")
    
    rasters = [rasterio_open(f) for f in files]
    
    # Merge the rasters into a single raster
    mosaic, out_trans = merge(rasters, nodata=0)  # Set nodata value within the uint16 range
    
    # Set output file path
    output_path = join(output_folder, f"new_mosaic_by_date_{date}.tif")
    
    # Write the mosaic to a GeoTIFF file
    with rasterio_open(output_path, 'w', driver='GTiff', height=mosaic.shape[1], width=mosaic.shape[2], count=mosaic.shape[0], dtype=mosaic.dtype, crs=rasters[0].crs, transform=out_trans, nodata=0) as dest:
        dest.write(mosaic)
    
    # Close the raster files
    for raster in rasters:
        raster.close()

    print(f"Mosaic {idx} for date {date} created at: {output_path}")
    print("=============================================")
