Meshing for a compound flood simulation can be very challenging.
Specifically, the small rivers in the watershed need to be resolved accurately and economically (in terms of mesh size) to ensure hydraulic connectivity,
and flood risk reduction structures such as levees and dams also need to be represented.

Here we dicuss the mesh generation for [STOFS3D](https://nauticalcharts.noaa.gov/updates/introducing-the-inland-coastal-flooding-operational-guidance-system-icogs/) v6.
A sample STOFS3D mesh is shown below (zoomed-in view on South Carolina):

![mesh stofs3d v6](../../assets/mesh-stofs-v6.png)

And here is an illustration of the evolution of STOFS3D meshes:

![mesh stofs3d v6](../../assets/mesh-stofs-versions.png)

The mesh generation is done in three steps:

1. Extracting thalwegs from DEM tiles.

2. Generate an SMS map that contains river arcs and intersection points.

3. Meshing in SMS.
