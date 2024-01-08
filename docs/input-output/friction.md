Bottom friction can be provided in three types of file - `drag.gr3`, or `rough.gr3` or `manning.gr3`. The ‘depth’ value in Gr3 file means $C_d$, bottom roughness in meters, or Manning’s $n$, respectively. The 3 files correspond to `nchi=0,1,-1`.

Bottom friction is a critical parameter in shallow area. Note that the bottom friction parameterizations are 
very different between 2D and 3D model, and so you cannot use same input. For details please read 
[this article](http://ccrm.vims.edu/yinglong/wiki_files/Report-ChezyFlow-Sept2011.pdf).

Further information can be found [here](../getting-started/grid-generation.md).
