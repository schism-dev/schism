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

## Advanced Parameterization
In the "Serial" and "Parallel" examples above, you may have noticed there are 3 mandatory Inputs:

| parameter | explanation |
| ----------- | ----------- |
| tif_fnames | a list of TIF file names. These TIFs should cover the area of interest and be arranged by priority (higher priority ones in front) |
| thalweg_shp_fname | name of a polyline shapefile containing the thalwegs |
| output_dir | must specify one. |

In addition to the mandatory inpouts, RiverMapper provides a few parameters to fine tune the output polylines or generate special features like levees or pseudo-channels.

| parameter | explanation |
| ----------- | ----------- |
| selected_thalweg | indices of selected thalwegs for which the river arcs will be sought. |
| output_prefix | a prefix of the output files, mainly used by the caller of this script; can be empty |
| mpi_print_prefix | a prefix string to identify the calling mpi processe in the output messages; can be empty |
| MapUnit2METER = 1 |  to be replaced by projection code, e.g., epsg: 4326, esri: 120008, etc. |
| river_threshold |  minimum and maximum river widths (in meters) to be resolved |
| outer_arc_positions | relative position of outer arcs, e.g., (0.1, 0.2) will add 2 outer arcs on each side of the river (4 in total), 0.1 \* riverwidth and 0.2 \* riverwidth from the banks. |
| length_width_ratio | a ratio of element length in the along-channel direction to river width; when a river is narrower than the lower limit, the bank will be nudged (see next parameter) to widen the river |
| i_close_poly | whether to add cross-channel arcs to enclose river arcs into a polygon |
| i_blast_intersection | whether to replace intersecting arcs (often noisy) at river intersections with scatter points (cleaner) |
| blast_radius_scale |  coef controlling the blast radius at intersections, a larger number leads to more intersection features being deleted |
| bomb_radius_coef |  coef controlling the spacing among intersection joints, a larger number leads to sparser intersection joints |
| snap_point_reso_ratio |  scaling the threshold of the point snapping; a negtive number means absolute distance value |
| snap_arc_reso_ratio |  scaling the threshold of the arc snapping; a negtive number means absolute distance value |
| i_DEM_cache  | Whether or not to read DEM info from cache.  Reading from original \*.tif files can be slow, so the default option is True |
| i_OCSMesh | Whether or not to generate outputs to be used as inputs to OCSMesh. |
| i_DiagnosticsOutput | whether to output diagnostic information |
| i_pseudo_channel | 0:  default, no pseudo channel, nrow_pseudo_channel and pseudo_channel_width are ignored; 1: fixed-width channel with nrow elements in the cross-channel direction, it can also be used to generate a fixed-width levee for a given levee centerline =2: implement a pseudo channel when the river is poorly defined in DEM
| pseudo_channel_width |  width of the pseudo channel (in meters) |
| nrow_pseudo_channel |  number of rows of elements in the cross-channel direction in the pseudo channel |

If their values are not explicitly specified, default values will be used.
You can also change the values of these parameters so that the output river map better fits your application.
For example, if you want to add a two pairs of outer arcs outside the river channel, you can do:
```python
    make_river_map(
        tif_fnames = ['./Inputs/DEMs/GA_dem_merged_ll.tif'],
        thalweg_shp_fname = './Inputs/Shapefiles/GA_local.shp',
        output_dir = './Outputs/',
        outer_arc_positions = (0.1, 0.2),
    )
```
where the last argument specifies the relative distance of the outer arcs to the main channel.

However, it may not be easy at first because there are many things to tune.
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

You can also peruse the "class methods" inside "Config_make_river_map" at the beginning of "make_river_map.py" to get a hang of what each parameter does.
Each preset (class method, or factory method) comes with a short description, for example,
"LooselyFollowRivers" means that "Small-scale river curvatures may not be exactly followed, but channel connectivity is still preserved".
To avoid duplication, the descriptions of the presets (which may be updated from time to time) are not listed in this manual.

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

Compare with the previous parallel sample, the above example adds a few lines of code.
It first imported the "Config_make_river_map", then created a default configuration,
and finally passed it to the parallel driver.
Since the default configuration is used, this example is essentially the same as the previous parallel example.
But if you don't want to start from the default configuration but from a preset,
and if you want to further modify the parameters based on the preset, you can do:
```python
    river_map_config = Config_make_river_map.LooselyFollowRivers()
    river_map_config.optional['i_DiagnosticOutput'] = True
    river_map_config.optional['i_real_clean'] = True
```
, without changing other parts of the code.

## Experimental Features
### Global arc cleaning
There is an optional input in Config_make_river_map called "i_realclean",
which defaults to False.
When turned on, an extra cleaning process will be applied on all arcs.
However, this process is not efficient enough for large applications,
specifically it increases the total running time 7 to 8 times.
Another issue is that it does not guarantee all river arcs are closed to form polygons,
which affects OCSMesh.
Updates in the near future will address the efficiency issue first.
