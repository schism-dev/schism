Some pre-proc and input files for the NOAA Pacific project

(1) DEM tiles: /sciclone/home20/whuang07/schism20/Pacific/DEM/RUNdir/DEM/dem_*.asc
     order in: /sciclone/home20/whuang07/schism20/Pacific/DEM/RUNdir/DEM/symlink_dems.pl 
(2) Interp script is interpolate_depth_structured2_mpi.f90
      vdatum(:) is used to convert vertical datum from local to MSL
(3) ~147 tide gauges are in: Pacific_tide_stations.bp
(4) new*.map are different versions of SMS maps
      new6: new4a-5 with finer resolution (8-10km) in ocean (@h=4km), 6km along US-Canada west coast. 2.4M nodes
      new7: from Wei. Refined west coast & AK from new6. ~2.5Mil nodes
      new10: new7 with a few arcs near Jap Trench removed; Seto Inland Sea, Guam edited to better acommodate
             tide gauges.
      new11: based on new10 with stations along west coast with bad performance beding refined based on the results from r10k
