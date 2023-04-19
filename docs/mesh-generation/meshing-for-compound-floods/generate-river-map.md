## Scripts
[RiverMapper](https://github.com/schism-dev/schism/tree/master/src/Utility/Grid_Scripts/Compound_flooding/RiverMapper) is available from the Utility folder in the SCHISM Git Repo.


## Usage
RiverMapper requires two inputs:

1. \*.tif, which are DEM tiles in lon/lat from one or more sources (such as CuDEM).

2. a shapefile with a 'LineString' type, which contains a 1D river network.

!!!note 
    The 1D river network can be any reasonable approximation of the thalwegs. It can be:
    <ul>
        <li>extracted from DEM using the method presented in the [previous section](./extract-thalweg.md)</li>
        <li>or duplicated from the river network of a hydrological model such as the National Water Model</li>
        <li>or manually drawn for [quick local touch-ups]()</li>
    </ul>

The outputs include:

- "total_river_arcs.map", which contains river arcs to be used for the final meshing in SMS.

- other \*.map files for diagnostic purposes.

A sample output looks like this:

![Sample river map](../../assets/sample-river-map.png) 


## Sample applications
To test the RiverMapper tool, you can begin by extracting the "RiverMapper_Samples/" directory from [RiverMapper_Samples.tar](http://ccrm.vims.edu/yinglong/feiye/Public/RiverMapper_Samples.tar).
It contains two subdirectories: "Serial" and "Parallel", providing sample applications for a smaller domain and a larger domain respectively.
Each subdirectory contains a sample Python script and the necessary input files. 

### Serial mode
For a small domain (covering one or two states), a direct function call to the serial "make_river_map" suffices.
See the sample script:
```
RiverMapper_Samples/Serial/sample_serial.py
```
, which reads:
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

Under the  directory "RiverMapper_Samples/Serial/", execute the serial script like this:
```
./sample_serial.py
```

### Parallel mode
For a large domain such as [STOFS3D Atlantic](https://nauticalcharts.noaa.gov/updates/introducing-the-inland-coastal-flooding-operational-guidance-system-icogs/),
a parallel driver is provided to automatically group thalwegs based on their parent tiles, then distribute the groups to parallel processors.

The sample parallel script is:
```
RiverMapper_Samples/Parallel/sample_parallel.py
```
, which reads:
```python
from mpi4py import MPI
import os
from RiverMapper.river_map_mpi_driver import river_map_mpi_driver

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

Under the directory "RiverMapper_Samples/Parallel/", execute the parallel script like this:

```
mpirun -n 20 ./sample_parallel.py
```



