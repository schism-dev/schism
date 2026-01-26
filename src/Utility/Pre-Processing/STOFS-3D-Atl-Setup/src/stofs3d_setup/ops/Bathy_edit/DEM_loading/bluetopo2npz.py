from pylib import convert_dem_format
import rasterio
from pathlib import Path
from glob import glob

def reproject_tif(input_geotiff, output_geotiff, dst_crs = 'EPSG:4326'):
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    # Open the source GeoTIFF file
    with rasterio.open(input_geotiff) as src:
        # Read the source CRS and calculate the transformation to WGS 84
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        
        # Update metadata for the output file
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        # Reproject and write to the output file
        with rasterio.open(output_geotiff, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)

def extract_band_tif(input_geotiff, output_geotiff, band_list=[]):
    if band_list == []:
        raise ValueError('band_list cannot be empty')

    with rasterio.open(input_geotiff) as src:
        # Read only the first band
        first_band = src.read(1)

        # Update the metadata for saving the single band
        meta = src.meta
        meta.update(count=1)

        # Write the first band to a new TIFF file
        with rasterio.open(output_geotiff, 'w', **meta) as dst:
            dst.write(first_band, 1)

if __name__ == "__main__":
    input_folder = '/sciclone/schism10/Hgrid_projects/DEMs/BlueTopo2/All_UTM_ZONES/'
    output_folder = '/sciclone/schism10/Hgrid_projects/DEMs/BlueTopo2/All_UTM_ZONES_lonlat_npz/'

    tif_files = [Path(file) for file in glob(f'{input_folder}/BlueTopo*.tiff')]
    for i, tif_file in enumerate(tif_files):
        band1_tif = f"{tif_file.parent}/{tif_file.stem}.band1.tif"
        extract_band_tif(tif_file, band1_tif, band_list=[1])

        projected_geotiff = f"{output_folder}/{tif_file.stem}.ll.tif"
        reproject_tif(band1_tif, projected_geotiff, dst_crs='EPSG:4326')

        output_fname = f"{output_folder}/{tif_file.stem.replace('.', '_')}_{i+1}.npz"
        S = convert_dem_format(projected_geotiff, sname=output_fname)

    pass