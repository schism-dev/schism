# Convert Datum from NAVD88 to xGEOID20b:

## 1. gen.py: prepare ascii file for NOAA's vertical Datum transformation software
    inputs:
      - chea_del_bay.reg 
      - stofs3d_inland2.reg 
      - hgrid.ll (after loading bathymetry, setting levee/feeder depth)

    outputs:
      - hgrid_nochesdel_navd_positiveup.txt
      - hgrid_chesdel_navd_positiveup.txt

## 2. VDatum:
    2.1 download the latest full vdatum package from https://vdatum.noaa.gov/docs/datums.html
        As of 02/24/2023, the version is vdatum_all_20221116.zip
    2.2 unzip the package and enter into the vdatum_all_20221116/vdatum/ folder, submit job.sh

    Outputs will be saved to a new created folder result/, and rename the file as:
        hgrid_stofs3d_inland_ches_del_New.txt
        hgrid_stofs3d_inland_1_New.txt
        hgrid_stofs3d_inland_2_New.txt
        hgrid_stofs3d_inland_3_New.txt
        hgrid_stofs3d_inland_4_New.txt
    
## 3. replace_depth.py
    inputs:
        hgrid_stofs3d_inland_ches_del_New.txt
        hgrid_stofs3d_inland_1_New.txt
        hgrid_stofs3d_inland_2_New.txt
        hgrid_stofs3d_inland_3_New.txt
        hgrid_stofs3d_inland_4_New.txt
        - hgrid.ll (same file as in step 1)

    outputs:
         hgrid_xGEOID20b.gr3
