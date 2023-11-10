# **Overview**
Here we use the SECOFS domain to illustrate the meshing procedure of a SMS project that has both manual polygons and automatically generated watershed river polygons.

<a href="../../assets/SECOFS_domain.jpg">
<img src="../../assets/SECOFS_domain.jpg" alt="drawing" width="500"/>
</a>

The key is to organize arcs and polygons in different SMS map coverages and merge them properly to avoid unwanted intersections.


# **Preparation**
The SECOFS SMS project has the following components:

![SMS project](../../assets/sms-project.jpg) 

For a complex meshing project, it is imperative to keep different features in seperate maps.
Merging individual maps should be done as the last step.

In this project, the purpose of each SMS map coverage is as follows:

- lbnd: land boundary roughly along the 10-m contour.

- levee\*: Levees from National Levee Database.

- auto\*arcs: Automatically generated river arcs, clipped inside the watershed region (see [clipping the river
  arcs](https://vims0-my.sharepoint.com/personal/feiye_vims_edu/_layouts/OneNote.aspx?id=%2Fpersonal%2Ffeiye_vims_edu%2FDocuments%2FNotebooks%2FNOAA%20TWL%20project%E2%80%8B&wd=target%28Highlights.one%7C9F3230A1-12F7-4F9A-8890-CE493995252A%2FSMS%20notes%7C95E15132-42B8-4581-A263-4F4576C9A5FE%2F%29onenote:https://vims0-my.sharepoint.com/personal/feiye_vims_edu/Documents/Notebooks/NOAA%20TWL%20project​/Highlights.one#SMS%20notes&section-id={9F3230A1-12F7-4F9A-8890-CE493995252A}&page-id={95E15132-42B8-4581-A263-4F4576C9A5FE}&end))

- coast: This map includes all manually made polygons (including quad patches) and the coastline.

- merge coverge: Merged map of all individual maps above.


!!! Attention
    Never edit the merged map directly, unless you are willing to take the time to exactly reproduce the changes in all individual maps (which is more time consuming in most cases, so don't do it).

    For a complex project, breaking the consistency between the merged map and the individual maps makes it practically impossible for any significant changes in the future.


# **Clean the merged map**

Since SMS 13.3, the cleaning is very efficient even for a map with tens of thousands of arcs. However, you may need to clean the merged map multiple times until there is no more warning messages.
For STOFS3D or SECOFS, use with the following parameters: 

![Clean SMS map](../../assets/mesh-clean-map.jpg)


# **Build polygons**

![Build polygons](../../assets/mesh-build-polygons.png)

Make sure there is no warning message after the polygons are built; otherwise, go back to the cleaning step and try again.


# **Finalize the mesh**
Delete the disjoint nodes:

![disjoint nodes](../../assets/mesh-disjoint-nodes.png)

Renumber nodes:

![renumber nodes](../../assets/mesh-renumber-nodes.png)

Save the mesh as \*.2dm:

![renumber nodes](../../assets/mesh-save.png)


# ** Mesh quality and skew elements **
Although the script tries to optimize the convergence of river arcs at river intersections
, skew elements may occasionally occur.
![river-intersections](../../assets/mesh-intersections.jpg)

Normally these small 'bad' elements do not affect the efficiency or stability of SCHISM simulations. 
However, it's advisable to edit out extremely small/skew elements.
You may remedy them using the following script provided in [RiverMapper](https://github.com/schism-dev/RiverMeshTools/tree/main/RiverMapper):
```
RiverMeshTools/RiverMapper/RiverMapper/improve_hgrid.py
```
The script automatically checks for small and skew elements,
then improve the mesh by merging, removing, or relaxing the problematic elements.
As a general rule of thumb, SCHISM can comfortably handle elements >= $1m^2$ ($10^{-10}$ in lon/lat), and skewness<=60. 
You can specify these parameters in the script; see the "__main__" section of the script for usage.

Alternatively, you can use `ACE/xmgredit5` or SMS to check and fix mesh quality manually.

!!! Attention
    Avoid overly depending on mesh adjustments in the final stage, as the underlying issue often originates from the SMS coverages or the merging process.
    In these instances, it is best to address the problems directly within the SMS coverages.


# **Additional notes for STOFS/SECOFS developers**
These notes are primarily intended as a reference for developers, although other users might also find parts of the content beneficial.


## Strategy for editing the coastal coverage
Edit the coastal coverage when local refinements are needed in an area of interest.

Rules to follow (for consistency among STOFS3D/SECOFS developers):

- The manual arcs should at least be better than the auto arcs; use the shapefile of the auto arcs as a background in SMS to help you decide.
  In particular, this means you need to accommodate for tributaries when placing a manual polygon for a main channel, lake, or estuary:

<img src="../../assets/accommodate_trib.jpg" alt="drawing" width="700"/>

- No need to create buffers for manual polygons; the [clipping procedure](https://vims0-my.sharepoint.com/personal/feiye_vims_edu/_layouts/OneNote.aspx?id=%2Fpersonal%2Ffeiye_vims_edu%2FDocuments%2FNotebooks%2FNOAA%20TWL%20project%E2%80%8B&wd=target%28Highlights.one%7C9F3230A1-12F7-4F9A-8890-CE493995252A%2FSMS%20notes%7C95E15132-42B8-4581-A263-4F4576C9A5FE%2F%29onenote:https://vims0-my.sharepoint.com/personal/feiye_vims_edu/Documents/Notebooks/NOAA%20TWL%20project​/Highlights.one#SMS%20notes&section-id={9F3230A1-12F7-4F9A-8890-CE493995252A}&page-id={95E15132-42B8-4581-A263-4F4576C9A5FE}&end) automatically creates a 50-meter buffer zone around all manual polygons.
!!! Attention
    Arcs not included in a polygon do not have a buffer; this setting is intentional.
  
- Unless you have a clear plan, don't mix features from other coverages (e.g., levees and auto arcs) in the coastal coverage, which is reserved for manual editing.
Doing so means that you assume the responsiblity of keeping the manual edits up-to-date with external changes that may be more easily done in other coverages.
For example, the National Levee Database add/remove features from time to time; the requirements on the auto arcs may change or there may be improvements in the RiverMapper tool itself.
In addition, putting auto arcs in this coverage means that you will spend extra time trying to satisfy the first rule of explicitly accommodating for the connectivity with tributaries.


## Set watershed resolution
To avoid over-refinement, most watershed polygons should have a specified mesh resolution, which can be set via the "constant paving" or "scalar paving density" option in each polygon's attributes.

Selecting the watershed polygons may involve much labor because many small polygons are generated after the river arcs are merged into the final map.
Instead of manual selection, use the ["select\*" map coverage](#select) to select all polygons more efficiently.

Activate the "select\*" coverage:

![Activate\_select](../../assets/mesh-activate-select.png)

Select the big polygon, right click on it, and click "Select intersecting objects":

![Intersect\_polygons](../../assets/mesh-intersect-polygons.png)

Intersect it with "merge coverage" with the following parameters:

![Intersect-parameters](../../assets/mesh-intersect-parameters.png)

It will take about 10 minutes to do the selection for the STOFS3D domain.

- For STOFS3D, right click on the selected polygons and set "mesh type" in "polygon attributes" to scalar paving density: 

![Scalar paving](../../assets/mesh-scalar-paving.png)

Set scalar options using the [watershed_resolution](#watershed_resolution) scatter dataset from the [preparation](#preparation) step:

![Scalar options](../../assets/mesh-scalar-options.png)

In addition, set the polygon attribute of the "island" between Chesapeake Bay and Delaware Bay as "None":

![Polygon None](../../assets/mesh-none.png)

- For SECOFS, "Constant paving" with the following parameters is used:

![const_paving](../../assets/const_paving.jpg)




