## Scripts
The scripts are available from this [Git Repo]().

## Usage
The inputs include the following:

1. DEM tiles in lon/lat from one or more sources (such as CuDEM).

2. Thalwegs specified in a line type shapefile.

For a small domain (one or two states), a direct function call to the serial "make_river_map" suffices, which looks like:

```python
make_river_map(
    tif_fnames = ['/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_dem_merged_utm17N.tif'],
    thalweg_shp_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/v4.shp',
    output_dir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/'
)
```

For a large domain such as [STOFS3D Atlantic](https://nauticalcharts.noaa.gov/updates/introducing-the-inland-coastal-flooding-operational-guidance-system-icogs/),
a parallel driver is provided to automatically group thalwegs based on their parent tiles, then distribute the groups to parallel processors.
The driver can be called like this:

```python
river_map_mpi_driver(
    dems_json_file=dems_json_file,
    thalweg_shp_fname=thalweg_shp_fname,
    output_dir=output_dir,
)
```

and executed like this:

```
mpirun -n 20 ./river_map_mpi_driver
```

The first argument takes a \*.json file that specifies multiple sets of DEMs (arranged from high priority to low priority), for example:
```json
{
    "CuDEM": {
        "name": "CuDEM",
        "glob_pattern": "/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/DEMs/CuDEM/Lonlat/*.tif",
        "file_list": [],
        "boxes": []
    },
    "CRM": {
        "name": "CRM",
        "glob_pattern": "/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/DEMs/CRM/Lonlat/*.tif",
        "file_list": [],
        "boxes": []
    }
}
```

## Products
The outputs of the script include the following files:

- "total_river_arcs.map", which contains river arcs to be used for the final meshing in SMS.

- "total_intersection_joints.map", which contains points used for identifying intersection elements to be relaxed in the [meshing stage](meshing-in-SMS.md#relax-crowded-elements-at-river-intersections).

, which can be imported to SMS directly.

The "total_river_arcs.map" needs to be cleaned to add feature points at line intersections and remove duplicated points.

This can be done in SMS, but it can be slow.
For example, cleaning the map for STOFS3D Atlantic can take hours and sometimes the intersection points are not generted correctly.
We have reported the problems to the SMS team and they are improving the cleaning process.

In the meantime, it is recommended to use Qgis to do the cleaning, see details [here]().


## Other comments

### Thalwegs

### Cleaning the river arcs
