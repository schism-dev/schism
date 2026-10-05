"""Tests for the shared USGS OGC request layer."""

import unittest
from unittest.mock import patch

from stofs3d_setup.ops.USGS.ogc import download_usgs


class UsgsOgcDownloaderTests(unittest.TestCase):
    def test_station_download_uses_ogc_station_identifier(self):
        feature = {
            "geometry": {"coordinates": [-73.75, 42.65]},
            "properties": {
                "monitoring_location_id": "USGS-01358000",
                "monitoring_location_name": "Hudson River at Albany",
                "time_series_id": "series-1",
                "time": "2017-01-01T00:00:00Z",
                "value": "10.5",
                "unit_of_measure": "ft3/s",
            },
        }

        with patch.object(
            download_usgs,
            "iter_ogc_features",
            return_value=[feature],
        ) as request:
            result = download_usgs.download_station_parameter(
                "01358000",
                "00060",
                "2017-01-01",
                "2017-01-02",
                session=object(),
            )

        self.assertEqual(result.loc[0, "station_id"], "01358000")
        self.assertEqual(result.loc[0, "value"], 10.5)
        query = request.call_args.args[2]
        self.assertEqual(query["monitoring_location_id"], "USGS-01358000")
        self.assertEqual(query["parameter_code"], "00060")

    def test_primary_only_group_download_filters_other_series(self):
        primary = {
            "id": "series-primary",
            "properties": {
                "monitoring_location_id": "USGS-01358000",
                "computation_identifier": "instantaneous",
                "begin_utc": "2010-01-01T00:00:00Z",
                "end_utc": None,
                "primary": "true",
            },
        }
        secondary = {
            "id": "series-secondary",
            "properties": {
                "monitoring_location_id": "USGS-01358000",
                "computation_identifier": "instantaneous",
                "begin_utc": "2010-01-01T00:00:00Z",
                "end_utc": None,
                "primary": "false",
            },
        }
        observations = [
            {
                "geometry": {"coordinates": [-73.75, 42.65]},
                "properties": {
                    "monitoring_location_id": "USGS-01358000",
                    "time_series_id": series_id,
                    "time": f"2017-01-01T{hour:02}:00:00Z",
                    "value": value,
                    "unit_of_measure": "ft3/s",
                },
            }
            for series_id, hour, value in [
                ("series-primary", 0, "1.0"),
                ("series-secondary", 1, "99.0"),
            ]
        ]

        with patch.object(
            download_usgs,
            "iter_ogc_features",
            side_effect=[[primary, secondary], observations],
        ):
            result = download_usgs.download_station_group_parameter(
                ["01358000"],
                "00060",
                "2017-01-01",
                "2017-01-02",
                session=object(),
                primary_only=True,
            )

        self.assertEqual(result["time_series_id"].tolist(), ["series-primary"])
        self.assertEqual(result["value"].tolist(), [1.0])


if __name__ == "__main__":
    unittest.main()
