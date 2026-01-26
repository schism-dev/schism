import geopandas as gpd
import rasterio
from pathlib import Path


# read an geotiff with RGB bands
tif_file = '/sciclone/schism10/Hgrid_projects/DEMs/BlueTopo/Savannah.tif'

output_file = f'{Path(tif_file).parent}/{Path(tif_file).stem}.xyz.tif'

with rasterio.open(tif_file) as src:
    # Read the red channel (assuming it's the first channel)
    elevation = src.read(1)

    # Define the new dataset's profile (similar to the source but with one band)
    profile = src.profile
    profile.update(
        dtype=rasterio.float32,
        count=1,
        compress='lzw'
    )

    # Write the elevation data to a new file
    with rasterio.open(output_file, 'w', **profile) as dst:
        dst.write(elevation, 1)

    pass