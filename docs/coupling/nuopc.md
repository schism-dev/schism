The National Unified Operational Prediction Capability (NUOPC) is a consortium of Navy, NOAA, and Air Force modelers and their research partners. It aims to advance the weather prediction modeling systems used by meteorologists, mission planners, and decision makers. NUOPC partners are working toward a common model architecture - a standard way of building models - in order to make it easier to collaboratively build modeling systems. To this end, they have developed the NUOPC Layer that defines conventions and a set of generic components for building coupled models using the Earth System Modeling Framework (ESMF).

SCHISM can be such a model component within an NUOPC-coupled system.  A so-called "cap" wraps SCHISM and exposes it via the NUOPC Application Programming Interface (API).  The interfaces exposed through the API are (1) import of fields, (2) export of fields, and (3) control structure.

## Obtaining and building the cap

The [NUOPC](nuopc.html)(esmf.html) cap are jointly hosted in a separate repository on [https://github.com/schism-dev/schism-esmf](https://github.com/schism-dev/schism-esmf).  It requires that the SCHISM core is built and pointed to by the environment variable `$SCHISM_BUILD_DIR` 

```
export SCHISM_ESMF_DIR=/my/path/to/schism-esmf
git clone https://github.com/schism-dev/schism-esmf.git $SCHISM_ESMF_DIR
cd $SCHISM_ESMF_DIR
make install-nuopc
```
The combined SCHISM and SCHISM-NUOPC libraries, as well as the NUOPC-compliant `Makefile` snippet `schism.mk` will be available in the `./lib` subdirectory.  

## Applications

The SCHISM NUOPC cap is currently used in the development version of the NOAA Environmental Modeling System (NEMS) [CoastalApp](coastalapp.html).  Within that system, SCHISM is coupled to components for the atmosphere (data or WRF), waves (data or WaveWatch III), and inland water (National Water Model).

## Reporting bugs or requesting features

The ESMF cap is still in development.  Please report any errors or annoyances in our bug tracker on https://github.com/schism-dev/schism-esmf/issues.  Also, please request features there, such as additional import or export fields.








