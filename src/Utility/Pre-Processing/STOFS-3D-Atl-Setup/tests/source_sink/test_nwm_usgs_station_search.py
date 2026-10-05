"""Regression tests for shared NWM-to-USGS reach traversal."""

import unittest

import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.nwm_usgs_station_search import (
    Vsource,
    find_usgs_along_nwm,
)


class NwmUsgsStationSearchTests(unittest.TestCase):
    def test_search_finds_first_gage_in_requested_direction(self):
        network = pd.DataFrame(
            {
                "Length": [100.0, 100.0, 100.0],
                "gages": [None, "upstream", "downstream"],
                "from": [[2], [], []],
                "to": [4, 3, 0],
            },
            index=[3, 2, 4],
        )
        source = Vsource(
            hgrid_ie=17,
            nwm_fid=3,
            df=pd.DataFrame(
                {
                    "datetime": pd.to_datetime(["2017-01-01"]),
                    "Data": [1.0],
                }
            ),
        )

        find_usgs_along_nwm(
            iupstream=True,
            starting_seg_id=source.nwm_fid,
            vsource=source,
            nwm_shp=network,
        )
        find_usgs_along_nwm(
            iupstream=False,
            starting_seg_id=source.nwm_fid,
            vsource=source,
            nwm_shp=network,
        )

        self.assertEqual(
            [station.st_id for station in source.usgs_st],
            ["upstream", "downstream"],
        )
        self.assertEqual(
            [station.nearby_nwm_fid for station in source.usgs_st],
            [2, 4],
        )
        self.assertEqual(source.upstream_path._seg_id, [3, 2])
        self.assertEqual(source.downstream_path._seg_id, [3, 4])


if __name__ == "__main__":
    unittest.main()
