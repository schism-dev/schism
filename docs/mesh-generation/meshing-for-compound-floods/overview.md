Meshing for a compound flood simulation can be very challenging.
It is crucial to accurately and economically resolve the small rivers in the watershed to ensure hydraulic connectivity while minimizing the mesh size.
Additionally, flood risk reduction structures such as levees and dams also need to be represented in the mesh.


The mesh of [STOFS3D](https://nauticalcharts.noaa.gov/updates/introducing-the-inland-coastal-flooding-operational-guidance-system-icogs/) has gone through major upgrades from STOFS3D v4 to v6,
which serves as a good illustration of different meshing techniques:

![mesh stofs3d v6](../../assets/mesh-stofs-versions.png)

In the v4 mesh, small watershed rivers are resolved by a 1D representation, using 1D river segments from the National Water Model as feature arcs to guide the meshing.

In the v5 mesh, characteristic contours of the DEM are extracted and used to represent river banks and thalwegs.
However, extensive manual cleaning is required due to the potential messiness of these contours.

In the v6 mesh, river arcs are generated automatically with
[pyDEM](https://github.com/schism-dev/RiverMeshTools/tree/main/pyDEM) and
[RiverMapper](https://github.com/schism-dev/RiverMeshTools/tree/main/RiverMapper),
both maintained in the [RiverMeshTools repository](https://github.com/schism-dev/RiverMeshTools).
Here is another sample v6 mesh that is zoomed in on the complex river network in South Carolina:
![mesh stofs3d v6](../../assets/mesh-stofs-v6.png)

The meshing procedure used in v6 is preferred due to its minimal need for manual labor and its ability to produce the best mesh quality out of the three versions.
The DEM-based procedure includes three steps:

1. [Extracting thalwegs from DEM tiles](extract-thalweg.md)

2. [Generate an SMS map that contains river arcs](generate-river-map.md)

3. [Meshing in SMS](meshing-in-SMS.md)

An existing channel network, such as NHD flowlines, can replace DEM-derived thalwegs. See the [NHD workflow](special-case-utilizing-nhd.md) for both real-DEM and NHD-area raster inputs.

## Publication
[A parallel Python-based tool for meshing watershed rivers at continental scale](https://www.sciencedirect.com/science/article/abs/pii/S1364815223001172)
