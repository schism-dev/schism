MOSSCO, the "Modular System for Shelves and Coasts" is a framework for coupling
processes or domains that are originally developed in standalone numerical models.
The software MOSSCO implements this infrastructure in the form of a library of
components and couplers, and of example coupled applications. The components
"wrap" external models used in coastal and shelf sciences; these wrapped components are then coupled
to each other in the Earth System Modeling Framework (ESMF).

The [SCHISM ESMF](esmf.html) cap integrates with MOSSCO.

## Obtaining and building MOSSCO

```
export SCHISM_BUILD_DIR=/my/path/to/schism/build
export SCHISM_ESMF_DIR=/my/path/to/schism-esmf
export MOSSCO_DIR=/my/path/to/mossco

git clone https://git.code.sf.net/p/mossco/code $MOSSCO_DIR

cd $MOSSCO_DIR
make all install

bash ./build.sh -component "SCHISM"
```

## Using SCHISM as part of a MOSSCO coupled system

A simple preconfigured application is available in `$MOSSCO_DIR/examples/esmf/schism`. To build it, run `make` in that directory.  You can use the resulting executable as a drop-in replacement for SCHISM's standalone `pschism` executable, but you need to add the `mossco.cfg` resource file which overrides `param.nml` for control parameters of the coupled system, like start and stop time.

## Reporting bugs or requesting features

The integration of SCHISM into MOSSCO is still in development.  Please report any errors or annoyances in the bug tracker on [https://sourceforge.net/p/mossco/tickets/](https://sourceforge.net/p/mossco/tickets/). 

