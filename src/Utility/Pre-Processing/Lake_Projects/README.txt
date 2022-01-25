SMS maps and other info for lake projects.
1. Superior_Wei_v8.map
     Comments: this map is the most updated map with a focus on Duluth Harbor. It also resolves all the major harbors around Lake Superior and incorporated with new changes from Josh's mesh.
     Important maps included in this map file are:
     	1) Superior_v6_10m previous version
	2) Load_depth_1: used for generating a separate mesh to load bathymetry for barrier islands around Duluth Harbor (-2m) and some Nemadji river channel (0.5m), the gr3 file generated using this mesh is called /sciclone/home20/whuang07/schism10/GreatLakes/Grid/v9/Depth_regions_1.gr3
	3) Load_depth_2: used for generating a separate mesh to load bathymetry for Duluth Harbor, bathymetry is obtained from Navigational Chart, the gr3 file generated using this mesh is called /sciclone/home20/whuang07/schism10/GreatLakes/Grid/v9/Depth_regions_from_NVChart.ll
	The method for generating this map is to generate a buffered polygon using the existing feature arcs with 0m contour.
	Depth are scattered points manually specified based on the Navigational chart.
     	4) Superior_v6_10m_dup: excluding Duluth Harbor, a previous map used to generate Superior_v7_10m_0
	5) Superior_v7_10m_0: a map with focus on Duluth Harbor, without merging NWM streams
	6) NWM_streams_5: NWM streams with 4 cross channel elements, 5 parallel lines
	7) Harbors: several other harbors are resolved by this map
	8) Stream_buffer_60m: used for generating a separate mesh to edit bathymetry along NWM streams (depression is 2m)
	The gr3 file is /sciclone/home20/whuang07/schism10/GreatLakes/Grid/v9/Depth_regions_3_streams2m.ll
	9) Superior_v7_10m_with_NWM: final version of Superior_v7 (a previous version)
	10) Superior_v7_10m_0_dup: a duplicate of Superior_v7_10m_0 with changes, the changes is made to merge the new changes from the mesh by Josh/Chin
	11) Superior_from_Josh: obtained from Josh/Chin (original without changes)
	12 Superior_from_Josh_dup: a duplicate of Superior_from_Josh, prepared to be merged with Superior_v7_10m_0_dup
	13) Superior_v8_10m_0=Superior_from_Josh_dup+Superior_v7_10m_0_dup_Harbors, a mesh merged incorporating Josh's mesh and the previous version (v7), before merging NWM streams
        14) Superior_v8_10m_withNWM: Final map including NWM streams

2. FiveLakes_Wei_v1.map
3. Other files
        1) DuluthHarbor_scatter_depth_fromNVC.h5 
           DEM_from_Chart: a set of scatter points with depth obtained from NOAA's navigational chart

