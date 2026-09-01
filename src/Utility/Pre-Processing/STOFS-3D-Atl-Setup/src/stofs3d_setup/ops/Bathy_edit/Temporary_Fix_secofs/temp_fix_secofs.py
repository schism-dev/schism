import copy
import geopandas as gpd
import numpy as np
from pathlib import Path
from pylib import read

from stofs3d_setup.utils.projection import project_geodataframe, project_grid


def temp_fix_secofs(gd_ll, wdir, reference_hgrid_file):
    """
    Replace hgrid depths within the SECOFS domain using depths
    from a reference hgrid.

    Parameters
    ----------
    gd_ll : schism_grid
        Target hgrid in lon/lat coordinates. Modified in place and returned.
    wdir : str
        Directory containing secofs_domain_clipped_v1.shp.
    reference_hgrid_file : str
        Reference hgrid whose depths are used inside the SECOFS domain.
    """

    # Project target grid to the same projected CRS used for spatial selection
    gd_meters = copy.deepcopy(gd_ll)
    project_grid(gd_meters, 'epsg:4326', 'esri:102008')

    hg_points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(gd_meters.x, gd_meters.y),
        crs='esri:102008'
    )

    # Read SECOFS domain polygon
    secofs_domain_file = (
        Path(wdir).parent
        / 'DEM_loading'
        / 'secofs_domain_clipped_v1.shp'
    )

    secofs_gdf = gpd.read_file(secofs_domain_file)


    secofs_gdf = project_geodataframe(
        secofs_gdf,
        'esri:102008',
        'SECOFS domain polygon projection'
    )

    # Find hgrid nodes inside the SECOFS domain
    joined_gdf = gpd.sjoin(
        hg_points,
        secofs_gdf,
        how='inner',
        predicate='within'
    )

    idx = joined_gdf.index.unique().to_numpy()

    print(
        f'{len(idx)} hgrid nodes found inside the SECOFS domain.'
    )

    # Read reference hgrid
    reference_hgrid = read(reference_hgrid_file)

    # Make sure node coordinates and ordering are identical
    coordinates_match = (
        gd_ll.np == reference_hgrid.np
        and np.allclose(gd_ll.x, reference_hgrid.x, atol=1e-8, rtol=0.0)
        and np.allclose(gd_ll.y, reference_hgrid.y, atol=1e-8, rtol=0.0)
    )

    if not coordinates_match:
        raise ValueError(
            'Target and reference hgrids do not have identical '
            'node coordinates and ordering.'
        )

    # Replace depth
    print('Replacing depths inside the SECOFS domain using reference hgrid.')
    gd_ll.dp[idx] = reference_hgrid.dp[idx]

    return gd_ll
