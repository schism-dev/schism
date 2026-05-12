'''
Sample script to download adt using copernicusmarine.
Register at https://data.marine.copernicus.eu/register?redirect=%2Fproducts
and save login info using "copernicusmarine login"
'''


import copernicusmarine

year = 2014
if year <= 2022:
    dataset_id = "c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1D"  # 31 Dec 1992 to 6 Jun 2023
else:  # 2023, 2024, ...
    dataset_id = "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D"  # 31 Dec 2021 to present

# B.C.
copernicusmarine.subset(
    dataset_id = dataset_id,
    variables = ["adt", "sla"],
    start_datetime = f"{year-1}-12-01T00:00:00", end_datetime = f"{year+1}-01-03T23:59:59",
    minimum_longitude = -60, maximum_longitude = -53, minimum_latitude = 1, maximum_latitude = 55,  # ocean boundaries
)

# I.C.
copernicusmarine.subset(
    dataset_id = dataset_id,
    variables = ["adt", "sla"],
    start_datetime = f"{year-1}-12-01T00:00:00", end_datetime = f"{year-1}-12-01T23:59:59",
    minimum_longitude = -105, maximum_longitude = -50, minimum_latitude = 1, maximum_latitude = 55,  # whole domain
)