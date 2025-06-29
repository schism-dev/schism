The `sflux/` dir is required if `nws=2`. In this case, atmoserpic forcings include wind, air pressure and temperature, precipitation, humidity and longwave and shortwave fluxes. These are specified in the netcdf files inside `sflux/` dir, and conform to the NetCDF Climate and Forecast (CF) Metadata Convention 1.0.

There are 4 types of files in sflux/dir; see [this site](http://ccrm.vims.edu/yinglong/wiki_files/NARR/) for sample files.

1. sflux_inputs.txt: This asci file is a namelist
2. sflux_air_1.X.nc: netcdf files that have time (in days), wind speed at 10m above MSL (u,v), air temperature and pressure and specific humidity;
3. sflux_prc_1.X.nc: netcdf files that time (in days), have precipitation data;
4. sflux_rad_1.X.nc: netcdf files that have time (in days), downward long and short (solar) wave radiation fluxes.

Note that X denotes the time stack number (1,2,3,...). There is no limit on max.


!!!note "sflux_input.txt"
    sflux_input.txt has the following basic structure - 

    ```
    &sflux_inputs ! file name
    /
    ```
    All parameters inside this input are optional. Advanced users may consult the source code for a complete list of parameters. Additionally, see [sample_input](https://raw.githubusercontent.com/schism-dev/schism/master/sample_inputs/sflux_inputs.txt) for a detailed example file.

!!!note "sflux_air"
    ```
    netcdf sflux_air_1.10 {
    dimensions:
            nx_grid = 349 ;
            ny_grid = 277 ;
            time = UNLIMITED ; // (8 currently)
    variables:
            float time(time) ;
                    time:long_name = "Time" ;
                    time:standard_name = "time" ;
                    time:units = "days since 2001-01-01" ;
                    time:base_date = 2001, 1, 1, 0 ;
            float lon(ny_grid, nx_grid) ;
                    lon:long_name = "Longitude" ;
                    lon:standard_name = "longitude" ;
                    lon:units = "degrees_east" ;
            float lat(ny_grid, nx_grid) ;
                    lat:long_name = "Latitude" ;
                    lat:standard_name = "latitude" ;
                    lat:units = "degrees_north" ;
            float uwind(time, ny_grid, nx_grid) ;
                    uwind:long_name = "Surface Eastward Air Velocity (10m AGL)" ;
                    uwind:standard_name = "eastward_wind" ;
                    uwind:units = "m/s" ;
            float vwind(time, ny_grid, nx_grid) ;
                    vwind:long_name = "Surface Northward Air Velocity (10m AGL)" ;
                    vwind:standard_name = "northward_wind" ;
                    vwind:units = "m/s" ;
            float prmsl(time, ny_grid, nx_grid) ;
                    prmsl:long_name = "Pressure reduced to MSL" ;
                    prmsl:standard_name = "air_pressure_at_sea_level" ;
                    prmsl:units = "Pa" ;
            float stmp(time, ny_grid, nx_grid) ;
                    stmp:long_name = "Surface Air Temperature (2m AGL)" ;
                    stmp:standard_name = "air_temperature" ;
                    stmp:units = "K" ;
            float spfh(time, ny_grid, nx_grid) ;
                    spfh:long_name = "Surface Specific Humidity (2m AGL)" ;
                    spfh:standard_name = "specific_humidity" ;
                    spfh:units = "1" ;

    // global attributes:
                    :Conventions = "CF-1.0" ;
    }
    ```

!!!notes "sflux_prc"
    ```
    netcdf sflux_prc_1.10 {
    dimensions:
            nx_grid = 349 ;
            ny_grid = 277 ;
            time = UNLIMITED ; // (8 currently)
    variables:
            float time(time) ;
                    time:long_name = "Time" ;
                    time:standard_name = "time" ;
                    time:units = "days since 2001-01-01" ;
                    time:base_date = 2001, 1, 1, 0 ;
            float lon(ny_grid, nx_grid) ;
                    lon:long_name = "Longitude" ;
                    lon:standard_name = "longitude" ;
                    lon:units = "degrees_east" ;
            float lat(ny_grid, nx_grid) ;
                    lat:long_name = "Latitude" ;
                    lat:standard_name = "latitude" ;
                    lat:units = "degrees_north" ;
            float prate(time, ny_grid, nx_grid) ;
                    prate:long_name = "Surface Precipitation Rate" ;
                    prate:standard_name = "precipitation_flux" ;
                    prate:units = "kg/m^2/s" ;

    // global attributes:
                    :Conventions = "CF-1.0" ;
    }
    ```

!!!notes "sflux_rad"
    ```
    netcdf sflux_rad_1.10 {
    dimensions:
            nx_grid = 349 ;
            ny_grid = 277 ;
            time = UNLIMITED ; // (8 currently)
    variables:
            float time(time) ;
                    time:long_name = "Time" ;
                    time:standard_name = "time" ;
                    time:units = "days since 2001-01-01" ;
                    time:base_date = 2001, 1, 1, 0 ;
            float lon(ny_grid, nx_grid) ;
                    lon:long_name = "Longitude" ;
                    lon:standard_name = "longitude" ;
                    lon:units = "degrees_east" ;
            float lat(ny_grid, nx_grid) ;
                    lat:long_name = "Latitude" ;
                    lat:standard_name = "latitude" ;
                    lat:units = "degrees_north" ;
            float dlwrf(time, ny_grid, nx_grid) ;
                    dlwrf:long_name = "Downward Long Wave Radiation Flux" ;
                    dlwrf:standard_name = "surface_downwelling_longwave_flux_in_air" ;
                    dlwrf:units = "W/m^2" ;
            float dswrf(time, ny_grid, nx_grid) ;
                    dswrf:long_name = "Downward Short Wave Radiation Flux" ;
                    dswrf:standard_name = "surface_downwelling_shortwave_flux_in_air" ;
                    dswrf:units = "W/m^2" ;

    // global attributes:
                    :Conventions = "CF-1.0" ;
    }
    ```

Note that `sflux_rad` is only required if the heat exchange module is invoked via `ihconsv=1`, and `sflux_prc` is only required if the salt exchange module is invoked via `isconsv=1`.
Since a barotropic model cannot do heat/salt exchange properly, these two types of sflux inputs should not be used there. To impose rainfall in a  barotropic model,
 you may use the source/sink option `if_source` by converting rainfall rate into sources.


We have NARR sflux files from 1979-present, but cannot upload all of them to the web due to disk space limitation. You can find some samples at [http://ccrm.vims.edu/yinglong/wiki_files/NARR/](http://ccrm.vims.edu/yinglong/wiki_files/NARR/).

Two sources of data are allowed for each type of `.nc` files, and the relative priority is fixed by the file name. For instance `sflux_air_1.3.nc` might be blended with a file called `sflux_air_2.3.nc`. The ".3" component of the name represents the order of the file within the stack of provided input files. For instance, there might be a new file (`1`, `2`, `3`) produced every 12 hours in a forecast cycle.

## Interpolation and prioritization
Using air as an example, it is assumed that the file `sflux_air_2.1.nc` is more resolved or accurate than sflux_air_1.1.nc. The two will be blended in the model in a way that favors the ‘_2’ file. This blending of the fields is only adjustable in the code as described in notes below. The default in sflux_9c.F90 is a 99:1 blend of the ‘_2’ file to the ‘_1’ file.

As was remarked above, the files are arranged temporally in a stack of files starting with ".1". Given the sequence of forecasting and analysis, it is common for atmospheric files to overlap. A file might begin with a brief period of data assimilation plus a few days of forecast. SCHISM assumes that a new file indicates the injection of information, so when it encounters overlap, it advances to the later file.

## Using NARR files for your simulation (North America only)
First, make sure the NARR grid covers your `hgrid.ll` to ensure proper sptatial interpolation.

In your run directory, `mkdir sflux` and inside it, create symbolic links to the NARR files. e.g., if you run starts from June 10, 2004 and ends June 20, 2004, then

```
sflux_air_1.1.nc --> narr_air.2004_06_10.nc
sflux_air_1.2.nc --> narr_air.2004_06_11.nc
...
sflux_air_1.11.nc --> narr_air.2004_06_20.nc
sflux_air_1.12.nc --> narr_air.2004_06_21.nc # (extra day to account for time zone difference)
```

Similarly for `sflux_rad_*.nc` and `sflux_prc_*.nc`. As described above, the number "1" after "air_" denotes the first data set used, with the second set taking priority; you can use up to 2 sets in SCHISM (which combines them with some given weights set in sflux_subs.F90); we only use 1 set in this example.

## Global sflux - NCEP CFSR
This [web site](http://rda.ucar.edu/) has many global atmospheric model files, and a very useful source is the NCEP's CFSR. Note that you need to register yourself at the site first.

Simply select 'NCEP Climate Forecast System Reanalysis (CFSR)', and then CFSRv2. Under 'Data Access' you can get a sub-set (instead of global files). Then select all variable required by SCHISM:

- Wind u,v at 10m
- (Air) pressure reduced to MSL
- (Air) temperature at 2m AGL
- Specific humidity at 2m AGL
- Downward longwave radiation flux at ground level
- Downward shortwave radiation flux at ground level
- Precipitation rate at ground level

In the new window, select output format as netcdf, and remember to restrict the vertical heights for each variable. Note that pressure is given at a different grid than other variables and so you'd do it separately. You can select a region and computational grid and temporal resolution etc and then submit the request.

Once you have downloaded the (compressed) netcdf files, you can use these [3 matlab scripts](http://ccrm.vims.edu/yinglong/wiki_files/process_CFSR.tgz) to process them into `sflux_[air,rad,prc]*.nc`. Note that the files are bundled differently for difference period (e.g. before and after Oct. 1, 2011) and so you may need to modify those scripts slightly.

## Preparing your own sflux inputs
After familiarize yourself with the NARR files and their format, you may embark on creating your own nc files. The best way is to modify existing matlab scripts (`src/Utility/Sflux_nc/readnc*.m`) included in the source code bundle, which have extensively in-line comments to guide you along the way.

Since the time and space interpolation will be used to interpolate sflux info onto hgrid.ll at runtime, you need to make sure that the lon/lat grid in `sflux_*_1.*` covers `hgrid.ll`, and the union of time records in `sflux*.nc` covers the entire simulation period. The time zone info is given by `utc_start`, and you may pre-pend some records if this value is negative (eastern hemisphere) or append if it is positive. Even if `utc_start=0`, the code will need at least one time record beyond `rnday` for interpolation; simply duplicate the record at `rnday` if you do not have such info.

!!!notes "Wind convention"
    u-component is eastward, v-comp. is northward (normal math convention, not compass convention)

!!!notes "Additional files"
    - `windrot_geo2proj.gr3`: rotates winds in case they do not align with coordinate axes, i.e. lat/lon
    <!--- - `watertype.gr3`: 6 is clear water, 7 is muddiest. Search in `schism_init.F90` for different water types. Required only if ihconsv=1.-->

!!!notes "Some details from sflux_9c.F90"
    - Of all attributes in nc file, only 'base_date' is required. This is the time origin used in each file and the ‘time’ is then the offset (in days) from this origin. Note that the last number (hour) is NOT used and UTC is always assumed in each file. Use `utc_start` in `param.nml` to adjust time zone;
    - The grids for air, rad and prc can be different (but must be the same within each type and each source). Additional requirements for the structured grid in .nc: [lon,lat](nx,ny) give x,y coord., nx is # of pts in x. Suppose a node in the grid is given by (i,j) (1<=i<=nx), then **the quad (i,j), (i+1,j), (i+1,j+1,i,j+1) can be in either counter-clockwise or clockwise direction (but must be self consistent within each set - no mix-and-match)**;
    - Search for "relative_weight" (inside `netcdf_io`) to change relative weights of the 2 sources for air, rad and prc if needed. All weights must > 0!
    - in case of 2 sources/grids for a variable, use "1" as larger grid (i.e. encompassing hgrid.ll) and "2" as smaller grid. The code will calculate weights associated with the 2 grids, and if some nodes in hgrid.ll fall outside grid "2" the interpolation will be done on grid "1" only (see combine_sflux_data, in particular, bad_node_ based on area coordinates outside [0,1]). Both grids must start from stack 1 but may have different # of stacks for each variable. Within each nc file # of time steps can vary. The cumulative time window of '2' does not need to cover the entire simulation (code will use values from '1' only if '2' time is missing), but window of '1' must;
    - `air_1_max_window_hours` (etc) are set in netcdf_io to define the max time stamp (offset from time origin) within each nc file. Besides those in netcdf_io, max_file_times (max. #of time records in each nc file) in routine get_times_etc() may need to be adjusted as well.
