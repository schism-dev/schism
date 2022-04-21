The NOAA Environmental Modeling System (NEMS) Coastal Application "CoastalApp" is a NUOPC-based coupled system using the NEMS coupler.

The [SCHISM NUOPC](nuopc.html) cap is currently being integrated as an OCN component of the CoastalApp; this is available from a public repository [https://github.com/noaa-ocs-modeling/CoastalApp](https://github.com/noaa-ocs-modeling/CoastalApp)

## Obtaining and building CoastalApp

```
export COASTALAPP_DIR=/my/path/to/coastalapp
git clone https://github.com/noaa-ocs-modeling/CoastalApp $COASTALAPP_DIR
cd $COASTALAPP_DIR
git checkout feature/schism
git submodule update --init --recursive SCHISM NEMS 
bash ./build.sh -component "SCHISM"
```

You can add components like `WW3` or `ADCIRC`, and you may be required to choose a compiler or platform.  Consult `./build.sh -h` for help and further information.

## Reporting bugs or requesting features

The integration of SCHISM into CoastalApp is still in development.  Please report any errors or annoyances in our bug tracker on [https://github.com/schism-dev/CoastalApp/issues](https://github.com/schism-dev/CoastalApp/issues). 

