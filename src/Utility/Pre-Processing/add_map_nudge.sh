#nco must install with netcdf-4 support (e.g. vortex)
ncap2 -O  -s 'map_to_global_node[$node]=array(1,1,$node)' SAL_nu.nc SAL_nu.nc.new
ncap2 -O  -s 'map_to_global_node[$node]=array(1,1,$node)' TEM_nu.nc TEM_nu.nc.new
#ncrename -d one,ntracers SAL_nu.nc.new
#ncrename -d one,ntracers TEM_nu.nc.new
