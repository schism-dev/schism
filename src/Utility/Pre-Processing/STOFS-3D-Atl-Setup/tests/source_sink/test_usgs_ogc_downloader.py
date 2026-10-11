"""Tests for the shared USGS OGC request layer."""

import unittest
from unittest.mock import Mock, patch

import pandas as pd

from stofs3d_setup.ops.USGS.ogc import download_usgs
from stofs3d_setup.ops.USGS.ogc import source_sink


class UsgsOgcDownloaderTests(unittest.TestCase):
    def test_feature_coordinates_skip_missing_or_invalid_geometry(self):
        self.assertEqual(
            download_usgs.feature_lon_lat({"geometry": None}),
            (None, None),
        )
        self.assertEqual(
            download_usgs.feature_lon_lat(
                {"geometry": {"coordinates": [None, 42.0]}}
            ),
            (None, None),
        )

    def test_source_sink_station_inventory_uses_location_metadata(self):
        locations = pd.DataFrame(
            {
                "station_id": ["USGS-01358000", "USGS-01359000"],
                "station_name": ["Hudson River at Albany", "Other station"],
                "state": ["NY", "NY"],
                "lon": [-73.75, -73.8],
                "lat": [42.65, 42.7],
            }
        )
        with patch.object(
            source_sink,
            "discover_monitoring_locations",
            return_value=locations,
        ) as inventory, patch.object(
            source_sink,
            "discover_station_ids_with_parameter",
            return_value={"01358000"},
        ) as screen:
            result = source_sink.get_usgs_stations_from_state(
                states=["NY"],
                parameter="00060",
                start="2017-01-01",
                end="2017-01-05",
            )

        inventory.assert_called_once_with(states=["NY"])
        screen.assert_called_once_with(
            ["USGS-01358000", "USGS-01359000"],
            parameter="00060",
            start="2017-01-01",
            end="2017-01-05",
        )
        self.assertEqual(result.loc[0, "site_no"], "01358000")
        self.assertEqual(result.loc[0, "station_nm"], "Hudson River at Albany")
        self.assertEqual(len(result), 1)

    def test_station_parameter_screen_batches_and_checks_period(self):
        in_period = {
            "properties": {
                "monitoring_location_id": "USGS-01358000",
                "computation_identifier": "Instantaneous",
                "begin_utc": "2010-01-01T00:00:00Z",
                "end_utc": "2020-01-01T00:00:00Z",
            }
        }
        wrong_period = {
            "properties": {
                "monitoring_location_id": "USGS-01359000",
                "computation_identifier": "Instantaneous",
                "begin_utc": "2000-01-01T00:00:00Z",
                "end_utc": "2010-01-01T00:00:00Z",
            }
        }
        with patch.object(
            download_usgs,
            "iter_ogc_features_cql2",
            side_effect=[[in_period, wrong_period], []],
        ) as request:
            matched = download_usgs.discover_station_ids_with_parameter(
                ["01358000", "01359000", "01360000"],
                parameter="00010",
                start="2017-01-01",
                end="2017-01-05",
                session=object(),
                batch_size=2,
            )

        self.assertEqual(matched, {"01358000"})
        self.assertEqual(request.call_count, 2)
        first_filter = request.call_args_list[0].args[2]
        self.assertEqual(
            first_filter["args"][0],
            {
                "op": "in",
                "args": [{"property": "parameter_code"}, ["00010"]],
            },
        )
        self.assertEqual(
            first_filter["args"][1]["args"][1],
            ["USGS-01358000", "USGS-01359000"],
        )

    def test_station_parameter_screen_propagates_metadata_errors(self):
        with patch.object(
            download_usgs,
            "iter_ogc_features_cql2",
            side_effect=RuntimeError("metadata query failed"),
        ):
            with self.assertRaisesRegex(RuntimeError, "metadata query failed"):
                download_usgs.discover_station_ids_with_parameter(
                    ["01358000"],
                    parameter="00060",
                    start="2017-01-01",
                    end="2017-01-05",
                    session=object(),
                )

    def test_monitoring_location_inventory_queries_by_state_code(self):
        feature = {
            "id": "USGS-01358000",
            "geometry": {"coordinates": [-73.75, 42.65]},
            "properties": {
                "agency_code": "USGS",
                "monitoring_location_number": "01358000",
                "monitoring_location_name": "Hudson River at Albany",
            },
        }
        feature_without_coordinates = {
            "id": "USGS-01358001",
            "geometry": None,
            "properties": {
                "agency_code": "USGS",
                "monitoring_location_number": "01358001",
            },
        }
        with patch.object(
            download_usgs,
            "iter_ogc_features",
            return_value=[feature_without_coordinates, feature],
        ) as request:
            result = download_usgs.discover_monitoring_locations(
                states=["NY"],
                session=object(),
            )

        self.assertEqual(result.loc[0, "station_id"], "01358000")
        self.assertEqual(result.loc[0, "station_name"], "Hudson River at Albany")
        self.assertEqual(result.loc[0, "lon"], -73.75)
        self.assertEqual(result.loc[0, "lat"], 42.65)
        self.assertEqual(request.call_args.args[1], "monitoring-locations")
        self.assertEqual(request.call_args.args[2], {"state_code": "36"})

    def test_cql2_inventory_query_preserves_filter_across_pages(self):
        filter_expression = {"op": "=", "args": [{"property": "x"}, 1]}
        with patch.object(
            download_usgs,
            "request_json",
            side_effect=[
                {
                    "features": [{"id": "first"}],
                    "links": [{"rel": "next", "href": "?offset=1"}],
                },
                {"features": [{"id": "second"}], "links": []},
            ],
        ) as request:
            result = download_usgs.iter_ogc_features_cql2(
                session=object(),
                collection="time-series-metadata",
                filter_expression=filter_expression,
                limit=10,
            )

        self.assertEqual([feature["id"] for feature in result], ["first", "second"])
        self.assertEqual(request.call_count, 2)
        self.assertEqual(
            request.call_args_list[0].kwargs["json_body"],
            filter_expression,
        )
        self.assertEqual(
            request.call_args_list[1].kwargs["json_body"],
            filter_expression,
        )

    def test_cql2_request_uses_usgs_filter_content_type(self):
        response = Mock()
        response.status_code = 200
        response.json.return_value = {"features": []}
        session = Mock()
        session.headers = {}
        session.post.return_value = response

        payload = download_usgs.request_json(
            session,
            "https://example.test/items",
            params={"f": "json"},
            json_body={"op": "=", "args": [{"property": "x"}, 1]},
        )

        self.assertEqual(payload, {"features": []})
        self.assertEqual(
            session.post.call_args.kwargs["headers"]["Content-Type"],
            "application/query-cql-json",
        )

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
