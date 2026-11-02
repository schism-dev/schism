Besides the following info, you may also search archived messages from the SCHISM mailing list on the same web where you registered yourself for the mailing list. Just log in, and select 'Archive' from the left of the screen, and then on right side of screen select 'Advanced search'. 

???faq "My results show platform/compiler/CPU dependency"
    Due to some intricate differences between different compilers/platforms some minor differences in results are expected. Bit-by-bit match of results using different numbers of CPUs is impossible also. But the differences should be small and should stabilize over time iteration.

???faq "Run crashed with a fatal.error error message "QUICKSEARCH: no intersecting edge....""
    First you need to viz the results before the crash (e.g. surface velocity) to see if anything is 
    outrageously wrong. Often times, errors such as this, related to the backtracking, come from
    other errors (e.g. NaN in inputs). 

    Very rarely, the error originates from the underflow issue in the backtracking/ELM (and the results otherwise 
    look reasonable). In the backtracking step, the code assumes that the path intersects sides of elements at 1 unique
    point that does not coincide with any node (cf. Fig. [1](./schism/eulerian-lagrangian-method.md#figure01)). Obviously this causes problem especially in some academic cases (e.g.
    the velocity aligns perfectly with grid line). To avoid this issue, the code nudges the starting position
    of backtracking (with parameter `btrack_nudge`) toward the starting element's centroid. Occasionally, this is still
    insufficient when the velocity is large on very fine meshes. One way out of this is to reduce `dtb_max` and `dtb_min`
    to get more accurate trajectory.

???faq "Run crashed with a fatal.error error message "bktrk_subs: overflow..." or "MAIN: nbtrk > mxnbt""
    The backtracking (ELM) in SCHISM is parallelized across MPI processes, and some arrays need to be allocated to store trajectories that exited the augmented domain. The dimension of those arrays is defined in mxnbt and mnbt, which are proportional to local # of side and # of levels, and the proportionality constants are defined in param.nml:

    `s1_mxnbt = 0.5`, `s2_mxnbt = 3.5`

    (Gradually) increasing these (to say 1, and 4) will solve the problem. Note that these constants only affect memory consumption and not accuracy.

???faq "How to do a tidal simulation with SCHISM?"
    The simplest way is to use SCHISM 2D. If you have to use SCHISM 3D, make sure the vertical grid is reasonable.

???faq "My run crashed; how can I find out why?"
    The best way to find out the reason for a crash is to visualize the surface velocity with VisIT. Usually you may see some large/noisy velocity somewhere, which may give you some hints on grid or forcing issues etc.
    
    Sometimes you want to visualize the problem right before the crash. You cannot use `autocombine_MPI_elfe.pl` as the last stack of output is incomplete. But you can use the FORTRAN combine script (e.g., `combine_output*`) to directly combine an incomplete stack. Just follow the instruction in the header of this script to prepare the inputs and run. Then visualize the combined outputs.

???faq "How to impose river discharge if the depth is negative there?"
    See [Grid generation](getting-started/grid-generation.md#grid-near-wetting-and-drying).

???faq "Run crashes with a dry boundary/depth error."
    SCHISM allows dry open-boundary nodes, but the entire open boundary segment cannot be dry at any given time if the velocity B.C. is imposed there. Either relocate the boundary or deepen some nodes to allow flow to come in. See [Grid generation](getting-started/grid-generation.md#grid-near-wetting-and-drying).

???faq "I have large velocity at open boundary"
    In hydrostatic models, the incoming velocity should be specified at open boundary. Over-specification i.e. elevation+velocity B.C. there usually avoids the problem, but the two B.C.’s must be largely consistent with each other. The velocity B.C. can be generated using global circulation model (HYCOM etc) plus global tidal model (FES etc).

???faq "My 3D velocity has noise or weird patterns"
    This can be caused by a number of reasons, but a common one is that when user switches from 2D to 3D, they forgot to set itur=3 to use turbulence closure.

???faq "fatal.error indicates “Could not find suitable element in input”"
    This is likely because your `sflux*.nc` is not prepared properly: either the grid in `.nc` does not cover `hgrid.ll`, or the structured grid in `.nc` has a mixture of clockwise and counter-clockwise orientations.

???faq "How to interpret elevation outputs from SCHISM?"
    Raw elevation outputs from SCHISM should be filtered in general to account for wetting and drying.  SCHISM's wetting and drying scheme is based on element (no partial element drying), and wet/dry flags at nodes/sides are derived from those at element. The dry flags can be found in outputs (e.g., "dryFlag*" in scribe IO). Since the elevation is defined at nodes, users should use "dryFlagNode" of the time record. If it is 1 (i.e. dry), then the elevation values should be ignored (as they are just placeholders for the next wet instance): values of 1000m or -1000m mean nothing in this case.

    Occasionally you may observe wet nodes having very large elevations. This may be due to the rare cases of rainfall falling on high ground, with no wet element around. So in addition to checking the wet/dry flag at a node, users should also check the wet/dry flags of its surrounding elements.  If surrounding  elements are all dry, then the node should be reclassified as dry, regardless of its original wet/dry flag. However, the large elevation values may be valid in other circumstances; e.g., if the local depth is -1000m (a mountain) and the elevation is 1000.1m, the total depth is only 10cm. If there is a wet element around this node, the elevation value is actually valid.

   In summary, follow these steps to interpret raw elevation outputs from SCHISM:

 - Check 'dryFlagNode'; if it is 1, then the node is dry at this time record and ignore the elevation
 - If dryFlagNode is 0 (wet), check the surrounding elements of this node. If all elements are dry at this time, then the node should be marked as dry at this time, and ignore the elevation value.
