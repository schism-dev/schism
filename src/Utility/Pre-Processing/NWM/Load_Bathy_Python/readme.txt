1). to convert a *asc file to *.npz file: use pylibs function convert_dem_fortmat 

2). use pload_depth.py to load bathymetry in parallel

   Input: 
     grd,grdout: names of orignal grid and grid with bathymetry loaded
     regions, rvalues: specify regions where you want to modified depth values 
     sdir: directory of DEM files (*.npz)
     qnode,nnode,ppn: (cluster name, # of nodes, # of core per node),choose cluster where you want to submit the script. 

