'''
Sample script to download adt using copernicusmarine.
Register at https://data.marine.copernicus.eu/register?redirect=%2Fproducts
and save login info using "copernicusmarine login"
'''


import copernicusmarine

copernicusmarine.subset(
  dataset_id = "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D",
  variables = ["adt", "sla"],
  start_datetime = "2024-04-01T00:00:00",
  end_datetime = "2024-04-11T23:59:59",
  minimum_longitude = -66,
  maximum_longitude = -60,
  minimum_latitude = 34,
  maximum_latitude = 38,
)
