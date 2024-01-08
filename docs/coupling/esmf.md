The Earth System Modeling Framework (ESMF) is a suite of software tools for developing high-performance, multi-component Earth science modeling applications. Such applications may include a few or dozens of components representing atmospheric, oceanic, terrestrial, or other physical domains, and their constituent processes (dynamical, chemical, biological, etc.). Often these components are developed by different groups independently, and must be “coupled” together using software that transfers and transforms data among the components in order to form functional simulations.

SCHISM can be such a model component within an ESMF-coupled system.  A so-called "cap" wraps SCHISM and exposes it via the ESMF Application Programming Interface (API).  The interfaces exposed through the API are (1) import of fields, (2) export of fields, and (3) control structure.

## Obtaining and building the cap

The [ESMF](esmf.html) and the [NUOPC](nuopc.html) caps are jointly hosted in a separate repository on [https://github.com/schism-dev/schism-esmf](https://github.com/schism-dev/schism-esmf).  It requires that the SCHISM core is built and pointed to by the environment variable `$SCHISM_BUILD_DIR` 

```
export SCHISM_ESMF_DIR=/my/path/to/schism-esmf
git clone https://github.com/schism-dev/schism-esmf.git $SCHISM_ESMF_DIR
cd $SCHISM_ESMF_DIR
make install-esmf
```
The combined SCHISM and SCHISM-ESMF libraries will be available in the `./lib` subdirectory.  

## Pre-configured executables

While the API exposed through the SCHISM-ESMF library is the core functionality of the cap, there are also several pre-configured simple coupled systems available for you to  `make` or use as a template: 

```
concurrent_esmf_test.F90  
multi_schism.F90
schism_driver_interfaces.F90
schism_pdaf.F90
triple_schism.F90
```

Have a look at the `examples` subdirectory!

## Applications

(1) The SCHISM ESMF cap is used to couple SCHISM to the Parallel Data Assimilation Framework leveraging the control structures of ESMF, see the example `schism_pdaf.F90`.

(2) The SCHISM ESMF cap is used in the Modular System for Shelves and Coasts, see  [MOSSCO](mossco.html).  Within that system, SCHISM can be flexibly coupled to components for the atmosphere, waves, ocean BGC, generic input/output, sediment. 

## Publications
Yu, H.-C., Zhang, Y. J., Nerger, L., Lemmen, C., Yu, J. C. S., Chou, T.-Y., Chu, C.-H., and Terng, C.-T.: Development of a flexible data assimilation method in a 3D unstructured-grid ocean model under Earth System Modeling Framework, EGUsphere [preprint], [https://doi.org/10.5194/egusphere-2022-114](https://doi.org/10.5194/egusphere-2022-114), 2022.

## Reporting bugs or requesting features

The ESMF cap is still in development.  Please report any errors or annoyances in our bug tracker on https://github.com/schism-dev/schism-esmf/issues.  Also, please request features there, such as additional import or export fields.








