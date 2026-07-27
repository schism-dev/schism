import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import argparse
from glob import glob
from pathlib import Path
from mpi4py import MPI
import time
import numpy as np
import geopandas as gpd
from netCDF4 import Dataset
import xarray as xr
from pylib import read


def pts_in_polygon(pts, polygon):
    hg_xy_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(pts[:, 0], pts[:, 1]), crs="EPSG:4326")
    idx = gpd.sjoin(hg_xy_gdf, polygon, how="inner", predicate="within").index.values
    return np.sort(idx).astype(np.int64)


def flood_fill_connected_wet_nodes(wet_mask, seed_mask, inp):
    """
    Perform a flood fill (BFS) to find all connected wet nodes starting from the seed nodes.
    """
    from collections import deque

    visited = np.zeros_like(wet_mask, dtype=bool)

    seeds = np.where(wet_mask & seed_mask)[0]
    q = deque(seeds)
    visited[seeds] = True

    while q:  # loop all seeds and all connected wet nodes
        i = q.popleft()  # get the next node index from the queue

        for nb in inp[i]:  # check all neighbors of the current node
            if nb < 0:  # skip padded -1 neighbors
                continue

            if wet_mask[nb] and not visited[nb]:
                q.append(nb)  # add unvisited wet neighbor to the queue
                visited[nb] = True

    return visited


def flood_fill_connected_wet_sides(wet_mask, seed_mask, isidenode):
    """
    Find wet nodes connected to wet seed nodes through actual mesh sides.

    A side is wet only when both of its endpoint nodes are wet.  Unlike
    ``flood_fill_connected_wet_nodes``, this does not connect nodes that share
    an element without sharing a side (for example, opposite nodes of a quad).

    Parameters
    ----------
    wet_mask : ndarray of bool, shape (nnode,)
        True for wet nodes.
    seed_mask : ndarray of bool, shape (nnode,)
        True for seed nodes.
    isidenode : ndarray of int, shape (nside, 2)
        Endpoint node indices for every mesh side.

    Returns
    -------
    ndarray of bool, shape (nnode,)
        True for wet nodes connected to at least one wet seed by wet sides.
    """
    from collections import deque

    visited = np.zeros_like(wet_mask, dtype=bool)

    seeds = np.where(wet_mask & seed_mask)[0]
    q = deque(seeds)
    visited[seeds] = True

    # Retain only sides whose two endpoint nodes are wet, then build the
    # adjacency needed by the BFS.  This deliberately excludes element
    # diagonals that may be present in node-neighbor connectivity.
    wet_sides = isidenode[np.all(wet_mask[isidenode], axis=1)]
    neighbors = [[] for _ in range(wet_mask.size)]
    for node1, node2 in wet_sides:
        neighbors[node1].append(node2)
        neighbors[node2].append(node1)

    while q:
        node = q.popleft()
        for neighbor in neighbors[node]:
            if not visited[neighbor]:
                visited[neighbor] = True
                q.append(neighbor)

    return visited


def flood_fill_connected_wet_elements(wet_mask, seed_mask, ine, elnode):
    """
    Find wet nodes connected to wet seeds through fully wet elements.

    An element is wet only when all of its valid nodes are wet.  The flood
    fill uses the node-to-element mapping ``ine`` to enter wet elements and
    then connects all nodes belonging to those elements.

    Parameters
    ----------
    wet_mask : ndarray of bool, shape (nnode,)
        True for wet nodes.
    seed_mask : ndarray of bool, shape (nnode,)
        True for seed nodes.
    ine : ndarray of int, shape (nnode, max_neighbor_elements)
        Neighboring element indices for each node, padded with negative values.
    elnode : ndarray of int, shape (nelement, max_nodes_per_element)
        Node indices for each element, padded with negative values.

    Returns
    -------
    ndarray of bool, shape (nnode,)
        True for wet nodes connected to a wet seed through fully wet elements.
    """
    from collections import deque

    visited_nodes = np.zeros_like(wet_mask, dtype=bool)
    visited_elements = np.zeros(elnode.shape[0], dtype=bool)

    seeds = np.where(wet_mask & seed_mask)[0]
    q = deque(seeds)
    visited_nodes[seeds] = True

    # Ignore negative padding when determining whether every real node of an
    # element is wet.  This supports both triangles and quadrilaterals.
    valid_element_nodes = elnode >= 0
    wet_elements = np.ones(elnode.shape[0], dtype=bool)
    for column in range(elnode.shape[1]):
        valid = valid_element_nodes[:, column]
        wet_elements[valid] &= wet_mask[elnode[valid, column]]

    while q:
        node = q.popleft()
        for element in ine[node]:
            if element < 0 or visited_elements[element] or not wet_elements[element]:
                continue

            visited_elements[element] = True
            nodes = elnode[element, valid_element_nodes[element]]
            new_nodes = nodes[~visited_nodes[nodes]]
            visited_nodes[new_nodes] = True
            q.extend(new_nodes)

    return visited_nodes


def get_seed_mask(hg, seed_region_file_list, output_dir, read_cache=False, write_cache=False):
    """
    Get a boolean mask for nodes that are within the seed region defined by the shapefiles.
    Also cache the seed mask to a NetCDF file for faster loading in future runs.

    Parameters
    ----------
    hg : hgrid
        The horizontal grid object.
    seed_region_file_list : list of str
        List of file paths to the shapefiles defining the seed regions.
    output_dir : str
        Directory where the cache file will be stored or read from.
    read_cache : bool, optional
        Whether to read the seed mask from a cache file if it exists. Default is False.
    """
    cache_success = False
    cache_file = f"{output_dir}/seed_mask_cache.nc"
    if read_cache:
        if os.path.exists(cache_file):
            ds_cache = xr.open_dataset(cache_file)
            seed_mask = ds_cache["seed_mask"].values.astype(bool)
            ds_cache.close()
            cache_success = True
            print(f"Loaded seed mask from cache file {cache_file}.")
        else:
            print(f"Warning: Cache file {cache_file} not found. Computing seed mask from shapefiles.")

    if not cache_success:
        seed_mask = np.zeros_like(hg.x, dtype=bool)
        for seed_region_file in seed_region_file_list:
            seed_region = gpd.read_file(seed_region_file).to_crs("EPSG:4326")
            nodes_in_seed_zone = pts_in_polygon(hg.xy, seed_region)
            seed_mask[nodes_in_seed_zone] = True

        if write_cache:
            ds_cache = xr.Dataset(
                {
                    "seed_mask": (("nSCHISM_hgrid_node",), seed_mask.astype(np.int8)),
                },
                coords={
                    "nSCHISM_hgrid_node": np.arange(hg.np),
                },
            )
            ds_cache["seed_mask"].attrs = {
                "long_name": "seed node mask for flood fill",
                "flag_values": np.array([0, 1], dtype=np.int8),
                "flag_meanings": "0: not a seed node; 1: seed node",
            }
            ds_cache.to_netcdf(cache_file, encoding={"seed_mask": {"dtype": "int8", "zlib": True, "complevel": 4}})
            ds_cache.close()
            print(f"seed mask from shapefiles saved to cache file {cache_file}.")

    return seed_mask


def append_isolated_ponds_mask_to_nc(nc_file, pond_mask):
    """
    Append isolated pond mask directly into an existing SCHISM NetCDF file.

    Parameters
    ----------
    nc_file : str
        Existing NetCDF file to modify in place.
    pond_mask : ndarray
        shape = (ntime, nSCHISM_hgrid_node)
        1 = isolated pond node
        0 = otherwise
    """

    varname = "isolatedPondNode"

    pond_mask = pond_mask.astype(np.int8)

    with Dataset(nc_file, "a") as ds:
        # check if the dimensions match the existing NetCDF file
        nt, nnode = pond_mask.shape
        if nt != len(ds.dimensions["time"]):
            raise ValueError(
                f"Time dimension of pond_mask ({nt}) != NetCDF file ({len(ds.dimensions['time'])}).")
        if nnode != len(ds.dimensions["nSCHISM_hgrid_node"]):
            raise ValueError(
                f"Node dimension of pond_mask ({nnode}) != NetCDF file ({len(ds.dimensions['nSCHISM_hgrid_node'])}).")

        if varname in ds.variables:
            var = ds.variables[varname]
        else:
            var = ds.createVariable(
                varname,
                "i1",
                ("time", "nSCHISM_hgrid_node"),
                zlib=True,
                complevel=4,
            )

            var.long_name = "isolated pond node mask"
            var.flag_values = np.array([0, 1], dtype=np.int8)
            var.flag_meanings = (
                "0: not an isolated pond node; "
                "1: isolated pond node"
            )

        var[:, :] = pond_mask


def output_isolated_ponds_mask(pond_mask, time, output_file):
    """
    Parameters
    ----------
    pond_mask : ndarray
        shape = (ntime, np)
        1 = isolated pond node
        0 = otherwise
    """

    nSCHISM_hgrid_node = pond_mask.shape[1]

    ds = xr.Dataset(
        {
            "isolatedPondNode": (
                ("time", "nSCHISM_hgrid_node"),
                pond_mask.astype(np.int8),
            ),
        },
        coords={
            "time": time,
            "nSCHISM_hgrid_node": np.arange(nSCHISM_hgrid_node),
        },
    )

    ds["isolatedPondNode"].attrs = {
        "long_name": "isolated pond node mask",
        "flag_values": np.array([0, 1], dtype=np.int8),
        "flag_meanings": "0: not an isolated pond node; 1: isolated pond node",
    }

    ds.to_netcdf(
        output_file,
        encoding={
            "isolatedPondNode": {
                "dtype": "int8",
                "zlib": True,
                "complevel": 4,
            }
        },
    )


if __name__ == "__main__":
    # usage:
    # Serial, all output files:
    # python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output
    #
    # MPI, all output files distributed across ranks:
    # mpiexec -n 10 python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output
    #
    # Single output file, useful for job arrays or launching independent processes:
    # python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output --id_seq_no 1
    #
    # the schism_input_dir should contain hgrid.gr3 and seed region shapefiles (e.g., coastal.shp, lakes.shp)
    # the schism_output_dir should contain SCHISM output files (e.g., *out2d_*.nc or out2d_{id_seq_no}.nc)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # input arguments
    argparser = argparse.ArgumentParser(
        description="Append isolated pond node masks to SCHISM out2d NetCDF files.",
        epilog=(
            "Examples:\n"
            "  Serial, all output files:\n"
            "    python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output\n\n"
            "  MPI, all output files distributed across ranks:\n"
            "    mpiexec -n 10 python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output\n\n"
            "  Single output file, useful for job arrays or independent processes:\n"
            "    python filter_isolated_ponds.py --schism_input_dir /path/to/input --schism_output_dir /path/to/output --id_seq_no 1"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    argparser.add_argument("--schism_input_dir", help="Path to the SCHISM input directory", required=True)
    argparser.add_argument("--schism_output_dir", help="Path to the SCHISM output directory", required=True)
    argparser.add_argument(
        "--id_seq_no",
        help=(
            "Optional output file sequence number string (e.g., '1' or '0001'). "
            "If omitted, all *out2d_*.nc files in schism_output_dir are processed."
        ),
        required=False,
        default=None,
    )
    argparser.add_argument(
        "--flood-fill-connectivity",
        choices=("nodes", "sides", "elements"),
        default="sides",
        help=(
            "Connectivity used by the flood fill: 'sides' uses only wet "
            "sides from hg.isidenode (default), 'elements' uses only elements "
            "whose nodes are all wet, and 'nodes' uses hg.inp."
        ),
    )
    args = argparser.parse_args()

    schism_input_dir = Path(args.schism_input_dir)
    schism_output_dir = Path(args.schism_output_dir)
    id_seq_no = args.id_seq_no

    # ---------------------------------------------------------------------
    # schism_input_dir = Path('/sciclone/schism10/feiye/STOFS3D-v7.2/Remove_ponds/SCHISM_inputs2/')
    # schism_output_dir = Path('/sciclone/schism10/feiye/STOFS3D-v7.2/Remove_ponds/SCHISM_outputs2/')
    # ---------------------------------------------------------------------

    for path_name, path in [
        ("schism_input_dir", schism_input_dir),
        ("schism_output_dir", schism_output_dir),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"--{path_name} does not exist: {path}")
        if not path.is_dir():
            raise NotADirectoryError(f"--{path_name} is not a directory: {path}")

    # read hgrid and seed regions

    hgrid_file = f"{schism_input_dir}/hgrid.gr3"
    seed_region_file_list = [  # see original files at /sciclone/schism10/feiye/STOFS3D-v7.2/Remove_ponds/SCHISM_inputs/
        f"{schism_input_dir}/coastal.shp",
        f"{schism_input_dir}/lakes.shp",
        # f"{schism_input_dir}/total_river_polys.shp"
    ]
    diag_outputs = False
    # ---------------end input parameters----------------

    time_start = time.time()

    if id_seq_no is None:
        schism_output_files = sorted(glob(f"{schism_output_dir}/*out2d_*.nc"))
        if len(schism_output_files) == 0:
            raise FileNotFoundError(f"No SCHISM output files found in {schism_output_dir}")
        print(f"Found {len(schism_output_files)} SCHISM output files in {schism_output_dir}")
    else:
        schism_output_files = sorted(glob(f"{schism_output_dir}/out2d_{id_seq_no}.nc"))
        if len(schism_output_files) == 0:
            raise FileNotFoundError(
                f"No SCHISM output files found matching pattern: out2d_{id_seq_no}.nc in {schism_output_dir}"
            )
        print(f"Found {len(schism_output_files)} SCHISM output file(s) matching sequence {id_seq_no} in {schism_output_dir}")

    print(f"Rank {rank}: {schism_output_files}\n")

    t0 = time.time()
    hg = read(hgrid_file)
    hg.compute_all()  # ensure hg.inp, hg.isidenode, and hg.ine are available for flood fill
    print(f"Rank {rank}: Read hgrid and computed connectivity in {time.time() - t0:.2f} seconds")

    if rank == 0:
        t0 = time.time()
        # any wet nodes connected (directly or indrectly through wet nodes) to wet seeds are not isolated ponds
        seed_mask = get_seed_mask(hg, seed_region_file_list, output_dir=schism_input_dir, read_cache=False)
        print(f"Rank {rank}: Computed seed mask in {time.time() - t0:.2f} seconds")
    else:
        seed_mask = None
    seed_mask = comm.bcast(seed_mask, root=0)

    # allocate schism_output_files to processes
    my_files = [f for i, f in enumerate(schism_output_files) if i % size == rank]

    for my_file in my_files:
        # read raw outputs
        ds = xr.open_dataset(my_file)
        dryFlagNode = ds["dryFlagNode"].values
        if diag_outputs:
            elevation = ds["elevation"].values
            elevation_fill_value = ds["elevation"].encoding.get("missing_value")
            # elevation_fill_value = np.nan
            elevation[dryFlagNode == 1] = elevation_fill_value
            elevation_clean = elevation.copy()
        t = ds["time"].values
        ds.close()

        # identify isolated ponds
        nt = dryFlagNode.shape[0]
        pond_mask_total = np.zeros_like(dryFlagNode, dtype=bool)
        for it in range(nt):
            time0 = time.time()
            wet_mask = dryFlagNode[it, :] == 0
            # hg.plot(value=wet_mask.astype(int), fmt=1);  import matplotlib.pyplot as plt; plt.show()
            if args.flood_fill_connectivity == "sides":
                connected_wet_mask = flood_fill_connected_wet_sides(
                    wet_mask=wet_mask, seed_mask=seed_mask, isidenode=hg.isidenode,
                )
            elif args.flood_fill_connectivity == "elements":
                connected_wet_mask = flood_fill_connected_wet_elements(
                    wet_mask=wet_mask, seed_mask=seed_mask, ine=hg.ine, elnode=hg.elnode,
                )
            else:
                connected_wet_mask = flood_fill_connected_wet_nodes(
                    wet_mask=wet_mask, seed_mask=seed_mask, inp=hg.inp,
                )
            print(f"Rank {rank}: {my_file.split('/')[-1]} it={it}: flood fill time = {time.time() - time0:.2f} seconds")

            pond_mask = wet_mask & ~connected_wet_mask
            pond_mask_total[it, :] = pond_mask

            if diag_outputs:
                # import matplotlib.pyplot as plt
                # hg.plot(value=wet_mask.astype(int), fmt=1)
                # plt.show()

                # hg.dp = pond_mask.astype(int) * 2 + connected_wet_mask.astype(int)
                # hg.grd2sms(f"{schism_output_dir}/{my_file.split('/')[-1].replace('.nc', f'_pond_mask_it{it}.2dm')}")

                elevation_clean[it, pond_mask] = elevation_fill_value
                hg.dp = elevation[it, :]
                hg.grd2sms(
                    f"{schism_output_dir}/{my_file.split('/')[-1].replace('.nc', f'_elevation_it{it}.2dm')}")
                hg.dp = elevation_clean[it, :]
                hg.grd2sms(
                    f"{schism_output_dir}/{my_file.split('/')[-1].replace('.nc', f'_elevation_clean_it{it}.2dm')}")

        # Save isolated ponds mask to a new NetCDF file
        # output_isolated_ponds_mask(
        #     pond_mask=pond_mask_total, time=t,
        #     output_file=f"{schism_output_dir}/{my_file.split('/')[-1].replace('.nc', '_isolated_pond_mask.nc')}"
        # )

        # Append the isolated ponds mask to the original NetCDF file (overwriting it)
        append_isolated_ponds_mask_to_nc(my_file, pond_mask_total)

    print(f"Rank {rank}: Total time for processing {len(my_files)} files: {time.time() - time_start:.2f} seconds")
