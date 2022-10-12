The NOAA Environmental Modeling System (NEMS) Coastal Application "CoastalApp" is a NUOPC-based coupled system using the NEMS coupler.

The [SCHISM NUOPC](nuopc.html) cap is currently being integrated as an ocean component of the CoastalApp, replacing the twodimensional ADCIRC model; CoastalApp is available from a public repository [https://github.com/noaa-ocs-modeling/CoastalApp](https://github.com/noaa-ocs-modeling/CoastalApp), and its integration with SCHISM occurs in the `develop` branch

## Obtaining and building CoastalApp

```
export COASTALAPP_DIR=/my/path/to/coastalapp
git clone https://github.com/noaa-ocs-modeling/CoastalApp -b develop $COASTALAPP_DIR 
cd $COASTALAPP_DIR
git submodule update --init --recursive SCHISM/schism SCHISM/schism-esmf NEMS 
bash ./build.sh -component "SCHISM"
```

You can add components like `WW3` or `ADCIRC`, and you may be required to choose a `-compiler`  or `-platform`, or set environment variables like `PARMETIS` or `ESMFMKFILE`.  Consult `./build.sh -h` for help and further information.

## Reporting bugs or requesting features

The integration of SCHISM into CoastalApp is still in development.  Please report any errors or annoyances in the upstream bug tracker on [https://github.com/noaa-ocs-modeling/CoastalApp/issues](https://github.com/noaa-ocs-modeling/CoastalApp/issues). 

