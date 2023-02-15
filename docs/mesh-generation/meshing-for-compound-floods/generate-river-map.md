## Scripts
The group of scripts for generating river map are available from a sub-folder ([RiverMapper](https://github.com/schism-dev/schism/tree/master/src/Utility/Grid_Scripts/Compound_flooding/RiverMapper)) in the SCHISM Git Repo.


## Usage
The inputs include:

1. DEM tiles in lon/lat from one or more sources (such as CuDEM).

2. Thalwegs specified in a line type shapefile.

The outputs include:

- "total_river_arcs.map", which contains river arcs to be used for the final meshing in SMS.

- other \*.map files for diagnostic purposes.

A sample product is like this:

![Sample river map](../../assets/sample-river-map.png) 


## Sample applications
Two sample applications are available in a tar archive, which can be downloaded [here](http://ccrm.vims.edu/yinglong/feiye/Public/RiverMapper_Samples.tar).

For a small domain (covering one or two states), a direct function call to the serial "make_river_map" suffices:

```python
from RiverMapper.make_river_map import make_river_map


if __name__ == "__main__":
    '''
    A sample serial application of RiverMapper
    '''
    make_river_map(
        tif_fnames = ['./Inputs/DEMs/GA_dem_merged_ll.tif'],
        thalweg_shp_fname = './Inputs/Shapefiles/GA_local.shp',
        output_dir = './Outputs/',
        cache_folder = './Cache/'
    )
```

You can try it out by extracting "RiverMapper_Samples/Serial/" from [RiverMapper_Samples.tar](http://ccrm.vims.edu/yinglong/feiye/Public/RiverMapper_Samples.tar),
where you can find the sample python script and necessary input files.

For a large domain such as [STOFS3D Atlantic](https://nauticalcharts.noaa.gov/updates/introducing-the-inland-coastal-flooding-operational-guidance-system-icogs/),
a parallel driver is provided to automatically group thalwegs based on their parent tiles, then distribute the groups to parallel processors.
The driver is provided in RiverMapper's folder in SCHISM's git repository:

```
[your_schism_dir]/schism/src/Utility/Grid_Scripts/Compound_flooding/RiverMapper/RiverMapper/river_map_mpi_driver.py
```

You only need to change the parameters in the main function at the bottom of the script:

```python
if __name__ == "__main__":

# ------------------------- sample input ---------------------------
dems_json_file = './Inputs/DEMs/dems.json'  # specifying files for all DEM tiles
thalweg_shp_fname='./Inputs/Shapefiles/LA_local.shp'
output_dir = './Outputs/' +  f'{os.path.basename(thalweg_shp_fname).split(".")[0]}_{size}cores/'
# ------------------------- end input section ---------------------------

river_map_mpi_driver(
    dems_json_file=dems_json_file,
    thalweg_shp_fname=thalweg_shp_fname,
    output_dir=output_dir,
)

```

The first argument takes a \*.json file that specifies multiple sets of DEMs (arranged from high priority to low priority), for example:
```json
{
    "CuDEM": {
        "name": "CuDEM",
        "glob_pattern": "./Inputs/DEMs/CuDEM/*.tif",
        "file_list": [],
        "boxes": []
    },
    "CRM": {
        "name": "CRM",
        "glob_pattern": "./Inputs/DEMs/CRM/*.tif",
        "file_list": [],
        "boxes": []
    }
}
```

The driver script can be executed like this:

```
mpirun -n 20 ./river_map_mpi_driver
```


A sample application utilizing the parallel driver (RiverMapper_Samples/Parallel/) can be extracted from [RiverMapper_Samples.tar](http://ccrm.vims.edu/yinglong/feiye/Public/RiverMapper_Samples.tar)

