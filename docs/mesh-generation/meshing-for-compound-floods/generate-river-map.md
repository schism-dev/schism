## Scripts
[RiverMapper](https://github.com/schism-dev/RiverMeshTools/tree/main/RiverMapper) is available from the [RiverMeshTools repository](https://github.com/schism-dev/RiverMeshTools) under schism-dev.


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

The exact mpi command may vary based on your system.


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
from RiverMapper.make_river_map import Config_make_river_map

if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    # ------------------------- sample input ---------------------------
    dems_json_file = '/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v16/Inputs/dems.json'  # specifying files for all DEM tiles
    thalweg_shp_fname='/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v16/Shapefiles/CUDEM_Merged_for_v16.shp'
    output_dir = './Outputs/' +  f'{os.path.basename(thalweg_shp_fname).split(".")[0]}_{comm.Get_size()}-core/'
    # ------------------------- end input section ---------------------------
    river_map_mpi_driver(
        dems_json_file=dems_json_file,
        thalweg_shp_fname=thalweg_shp_fname,
        output_dir=output_dir,
        comm=comm
    )
```

The first argument takes a \*.json file that specifies multiple sets of DEMs, for example:
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
There are a few things to be noted:

* In stead of specifying a file list for the DEM tiles (as in the serial example),
you can either specify a glob pattern or a list of file names;
if both are specified, the script will automatically pool all files together and take the unique files.

* No need to manually specify the bounding "boxes", just leave it empty.

* Different DEM products should be arranged from high priority to low priority.
In the above example, "CuDEM" precedes "CRM", meaning that if a point is within CuDEM's spatial coverage,
then its elevation value is solely based on CuDEM.

> The order of the dictionary key is preserved as long as you are using Python 3.7 or a later version.
> We may replace the regular dictionary to OrderedDict in the future.

Under the directory "RiverMapper_Samples/Parallel/", execute the parallel script like this:

```
mpirun -n 20 ./sample_parallel.py
```

## Advanced Parameterization

### Mandatory inputs

In the "Serial" and "Parallel" examples above, you may have noticed there are 3 mandatory Inputs:

| parameter | explanation |
| ----------- | ----------- |
| tif_fnames (or a \*.json file if there are many tiles) | a list of TIF file names. These TIFs should cover the area of interest and be arranged by priority (higher priority ones in front) |
| thalweg_shp_fname | name of a polyline shapefile containing the thalwegs |
| output_dir | must specify one. |

### Optional inputs

In addition to the mandatory inpouts, RiverMapper provides a few parameters to fine tune the output polylines or generate special features like levees or pseudo-channels.

| parameter | type | explanation |
| ----------- | ----------- | ----------- |
| selected_thalweg | numpy array | Indices of a subset of thalwegs for which the river arcs will be sought; mainly used by the parallel driver |
| output_prefix | string | a prefix of the output files, mainly used by the caller of this script; can be empty |
| mpi_print_prefix | string | a prefix string to identify the calling mpi processe in the output messages; can be empty |
| MapUnit2METER = 1 | float |  no need to change; to be replaced by projection code, e.g., epsg: 4326, esri: 120008, etc. |
| river_threshold | float | minimum and maximum river widths (in meters) to be resolved |
| outer_arc_positions | a list/tuple of floats | relative position of outer arcs, e.g., (0.1, 0.2) will add 2 outer arcs on each side of the river (4 in total), 0.1 \* riverwidth and 0.2 \* riverwidth from the banks. |
| length_width_ratio | float |  a ratio of element length in the along-channel direction to river width; when a river is narrower than the lower limit, the bank will be nudged (see next parameter) to widen the river |
| i_close_poly | bool | whether to add cross-channel arcs to enclose river arcs into a polygon |
| i_blast_intersection | bool | whether to replace intersecting arcs (often noisy) at river intersections with scatter points (cleaner) |
| blast_radius_scale | float | coefficient controlling the blast radius at intersections, a larger number leads to more intersection features being deleted |
| bomb_radius_coef | float | coefficient controlling the spacing among intersection joints, a larger number leads to sparser intersection joints |
| snap_point_reso_ratio | float | scaling the threshold of the point snapping; a negtive number means absolute distance value |
| snap_arc_reso_ratio | float | scaling the threshold of the arc snapping; a negtive number means absolute distance value |
| i_DEM_cache  | bool | Whether or not to read DEM info from cache.  Reading from original \*.tif files can be slow, so the default option is True |
| i_OCSMesh | bool | Whether or not to generate outputs to be used as inputs to OCSMesh. |
| i_DiagnosticsOutput | bool | whether to output diagnostic information |
| i_pseudo_channel | int | 0:  default, no pseudo channel, nrow_pseudo_channel and pseudo_channel_width are ignored; 1: fixed-width channel with nrow elements in the cross-channel direction, it can also be used to generate a fixed-width levee for a given levee centerline =2: implement a pseudo channel when the river is poorly defined in DEM
| pseudo_channel_width | float | width of the pseudo channel (in meters) |
| nrow_pseudo_channel |int| number of rows of elements in the cross-channel direction in the pseudo channel |

You can change the values of these parameters so that the output river map better fits your application (otherwise default values are used).
For example, if you want to add two pairs of outer arcs that flank the main river channel, you can do:
```python
    make_river_map(
        tif_fnames = ['./Inputs/DEMs/GA_dem_merged_ll.tif'],
        thalweg_shp_fname = './Inputs/Shapefiles/GA_local.shp',
        output_dir = './Outputs/',
        outer_arc_positions = (0.1, 0.2),
    )
```
, where the last argument specifies the relative distances of the outer arcs to the main channel.
In this case two outer arcs will be placed to the left of the channel at distances of "$0.1 \times$ channel width" and "$0.2 \times$ channel width" from the left bank;
the same goes for the right bank, so 4 outer arcs in total.

However, tuning the parameters may not be easy at first because there are many of them.
To simplify the parameter configuration, RiverMapper offers a class called "Config_make_river_map",
which provides commonly used parameter presets.
You can use these presets to directly configure the parameters.
For example,
```python
    from RiverMapper.make_river_map import Config_make_river_map

    river_map_config = Config_make_river_map.LooselyFollowRivers()
    make_river_map(
        tif_fnames = ['./Inputs/DEMs/GA_dem_merged_ll.tif'],
        thalweg_shp_fname = './Inputs/Shapefiles/GA_local.shp',
        output_dir = './Outputs/',
        **river_map_config.optional,
    )

```
, where "LooselyFollowRivers" is a parameter preset defined in the "Config_make_river_map" class,
and the operator "\*\*" unpacks the optional parameters from the configuration.

You can also peruse the "class methods" inside "Config_make_river_map" at the beginning of "make_river_map.py" to familiarize with different parameter settings.
Each preset (class method, or factory method) comes with a short description in the code, for example,
"LooselyFollowRivers" means that "Small-scale river curvatures may not be exactly followed, but channel connectivity is still preserved".
The presets may be regularly updated, so their descriptions are not repeated here to avoid duplication and inconsistency.

Another utility of the "Config_make_river_map" class is to facilitate the parameter transfer from the parallel driver to the core routine.
The stardard way of setting the parameters in the parallel mode is as follows:

```python
from mpi4py import MPI
import os
from RiverMapper.river_map_mpi_driver import river_map_mpi_driver
from RiverMapper.make_river_map import Config_make_river_map

if __name__ == "__main__":
    comm=MPI.COMM_WORLD
# ------------------------- sample input ---------------------------
    dems_json_file = '/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v16/Inputs/dems.json'  # specifying files for all DEM tiles
    thalweg_shp_fname='/sciclone/schism10/Hgrid_projects/STOFS3D-V6/v16/Shapefiles/CUDEM_Merged_for_v16.shp'
    output_dir = './Outputs/' +  f'{os.path.basename(thalweg_shp_fname).split(".")[0]}_{comm.Get_size()}-core/'
    river_map_config = Config_make_river_map()
# ------------------------- end input section ---------------------------
    river_map_mpi_driver(
        dems_json_file=dems_json_file,
        thalweg_shp_fname=thalweg_shp_fname,
        output_dir=output_dir,
        river_map_config=river_map_config,
        comm=comm
    )
```

Comparing with the parallel [example](#parallel-mode) in the previous section, the above example adds a few lines of code that does the following:

* importing the "Config_make_river_map"
* creating a default configuration,
* passing the configration to the parallel driver.

Since the default configuration is used, this example is essentially the same as the previous parallel example.
But it facilitates further parameter tweaking, e.g.:
```python
    river_map_config = Config_make_river_map.LooselyFollowRivers()
    river_map_config.optional['i_DiagnosticOutput'] = True
    river_map_config.optional['i_real_clean'] = True
```
, without changing other parts of the code.

## Experimental Features
### Global arc cleaning
There is an optional input in Config_make_river_map called "i_real_clean",
which defaults to False.
When turned on, an extra cleaning process will be applied to all arcs.
However, this process is not efficient enough for large applications at the time being,
specifically it increases the total run time by 7 to 8 times.
Another issue is that it does not guarantee all river arcs are closed to form polygons,
which affects OCSMesh.
Updates in the near future will address the efficiency issue first.
