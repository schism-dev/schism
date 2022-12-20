# Convert Datum from NAVD88 to xGEOID20b:

## 1. gen_points.py: prepare ascii file for NOAA's vertical Datum transformation software
    inputs:
      - chea_del_bay.reg (VDatum treats contiguous United States and Cheasapeake/Delaware Bay as two differenct regions, therefore we needs two seperate ascii files otherwise it will fail)
      - hgrid.ll.new_dem_name (from load_bathymetry.py)
      - hgrid.ll.new_dem_id (from load_bathymetry.py)
      - hgrid.ll.new (from load_bathymetry.py)

    outputs:
      - hgrid_nochesdel_navd_positiveup.txt
      - hgrid_chesdel_navd_positiveup.txt

## 2. In VDatum software:
   2.1 Choose Region: Contiguous United States
       Vertical Information:
           Source: NAVD 88, Target: xGEOID20b(IGS14)
       ASCII File Conversion
           File name: hgrid_nochesdel_navd_positiveup.txt
       Ellipsoidal Transformation Epoch Inputs:
           the reference date of the input position: 2010.0
           the reference date of the output position: 2020.0

       output filename hgrid_nochesdel_navd_positiveup_xGEOID20b.txt

   2.2 Choose Region: Cheasapeake/Delaware Bay
       Vertical Information:
           Source: NAVD 88, Target: xGEOID20b(IGS14)
       ASCII File Conversion
           File name: hgrid_chesdel_navd_positiveup.txt
       Ellipsoidal Transformation Epoch Inputs:
           the reference date of the input position: 2010.0
           the reference date of the output position: 2020.0

       output filename hgrid_chesdel_navd_positiveup_xGEOID20b.txt

## 3. replace_depth.py
    inputs:
        - hgrid_nochesdel_navd_positiveup_xGEOID20b.txt
        - hgrid_chesdel_navd_positiveup_xGEOID20b.txt
        - hgrid.gr3 (or hgrid.npz)

    outputs:
         hgrid_xGEOID20b.gr3
