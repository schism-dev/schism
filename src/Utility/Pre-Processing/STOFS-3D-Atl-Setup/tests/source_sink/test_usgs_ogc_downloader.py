"""Tests for the shared USGS OGC request layer."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import pandas as pd
import numpy as np

from stofs3d_setup.ops.USGS.ogc import download_usgs
from stofs3d_setup.ops.USGS.ogc import source_sink


class UsgsOgcDownloaderTests(unittest.TestCase):
    def test_boundary_buffer_is_independent_of_station_batch(self):
        hgrid = SimpleNamespace(
            x=np.array([0.0, 0.0]),
            y=np.array([60.0, 0.0]),
            bndinfo=SimpleNamespace(nb=1, ibn=[np.array([0, 1])]),
        )
        high_latitude_station = np.array([[1.7, 60.0]])
        both_stations = np.array([[1.7, 60.0], [0.0, 0.0]])

        self.assertTrue(
            download_usgs._within_boundary_buffer(
                high_latitude_station, hgrid, 100.0
            )[0]
        )
        np.testing.assert_array_equal(
            download_usgs._within_boundary_buffer(
                both_stations, hgrid, 100.0
            ),
            [True, True],
        )

    def test_partial_current_cache_fetches_only_uncached_stations(self):
        dates = pd.date_range("2017-01-01", "2017-01-02", tz="UTC")
        cached_record = source_sink.GenericObsData(
            station_info={"id": "01358000"},
            df=pd.DataFrame({"date": dates, "value": [1.0, 2.0]}),
        )
        with TemporaryDirectory() as directory:
            cache_path = Path(directory) / "flow.pq"
            cache_path.touch()
            Path(f"{cache_path}.manifest.json").write_text(
                json.dumps({
                    "version": 1,
                    "quality_version": source_sink.CACHE_QUALITY_VERSION,
                    "parameter": "00060",
                    "stations": ["01358000", "01359000"],
                    "start": "2017-01-01T00:00:00Z",
                    "end": "2017-01-02T00:00:00Z",
                }),
                encoding="utf-8",
            )
            with (
                patch.object(source_sink, "_read_cache", return_value=[cached_record]),
                patch.object(source_sink, "make_session", return_value=object()),
                patch.object(
                    source_sink,
                    "download_station_group_parameter",
                    return_value=pd.DataFrame(),
                ) as download,
                patch.object(source_sink, "_write_cache"),
                patch.object(source_sink, "_write_cache_manifest"),
            ):
                records = source_sink.download_stations(
                    param_id="00060",
                    station_ids=["01358000", "01360000"],
                    cache_fname=str(cache_path),
                    datelist=dates,
                    allow_adaptive_retry=False,
                )

        self.assertEqual([record.station_info["id"] for record in records], ["01358000"])
        self.assertEqual(download.call_args.kwargs["station_ids"], ["01360000"])

    def test_legacy_cache_without_quality_filter_is_refetched(self):
        dates = pd.date_range("2017-01-01", "2017-01-02", tz="UTC")
        with TemporaryDirectory() as directory:
            cache_path = Path(directory) / "flow.pq"
            cache_path.touch()
            Path(f"{cache_path}.manifest.json").write_text(
                json.dumps({
                    "api": "usgs-waterdata-ogc-v1",
                    "parameter": "00060",
                    "stations": ["01358000"],
                    "start": dates[0].isoformat(),
                    "end": dates[-1].isoformat(),
                }),
                encoding="utf-8",
            )
            with (
                patch.object(source_sink, "_read_cache") as read_cache,
                patch.object(source_sink, "make_session", return_value=object()),
                patch.object(
                    source_sink,
                    "download_station_group_parameter",
                    return_value=pd.DataFrame(),
                ) as download,
                patch.object(source_sink, "_write_cache"),
                patch.object(source_sink, "_write_cache_manifest"),
            ):
                source_sink.download_stations(
                    param_id="00060",
                    station_ids=["01358000"],
                    cache_fname=str(cache_path),
                    datelist=dates,
                    allow_adaptive_retry=False,
                )

        read_cache.assert_not_called()
        self.assertEqual(download.call_args.kwargs["station_ids"], ["01358000"])

    def test_cached_temperature_below_cutoff_fails_required_station(self):
        dates = pd.date_range("2017-01-01", "2017-01-02", tz="UTC")
        bad_record = source_sink.GenericObsData(
            station_info={"id": "02228000"},
            df=pd.DataFrame({"date": dates, "value": [-118.0, -119.0]}),
        )
        with TemporaryDirectory() as directory:
            cache_path = Path(directory) / "temperature.pq"
            cache_path.touch()
            Path(f"{cache_path}.manifest.json").write_text(
                json.dumps({
                    "api": "usgs-waterdata-ogc-v1",
                    "quality_version": source_sink.CACHE_QUALITY_VERSION,
                    "parameter": "00010",
                    "stations": ["02228000"],
                    "start": dates[0].isoformat(),
                    "end": dates[-1].isoformat(),
                }),
                encoding="utf-8",
            )
            with (
                patch.object(source_sink, "_read_cache", return_value=[bad_record]),
                patch.object(source_sink, "make_session", return_value=object()),
                patch.object(
                    source_sink,
                    "download_station_group_parameter",
                    return_value=pd.DataFrame(),
                ) as download,
                patch.object(source_sink, "_write_cache"),
                patch.object(source_sink, "_write_cache_manifest"),
            ):
                with self.assertRaisesRegex(
                    source_sink.RequiredStationDownloadError, "02228000"
                ):
                    source_sink.download_stations(
                        param_id="00010",
                        station_ids=["02228000"],
                        required_station_ids=["02228000"],
                        cache_fname=str(cache_path),
                        datelist=dates,
                        allow_adaptive_retry=False,
                    )

        self.assertEqual(download.call_args.kwargs["station_ids"], ["02228000"])

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
            "discover_station_ids_with_parameter_in_bbox",
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
            bbox=(-73.8, 42.65, -73.75, 42.7),
            parameter="00060",
            start="2017-01-01",
            end="2017-01-05",
        )
        self.assertEqual(result.loc[0, "site_no"], "01358000")
        self.assertEqual(result.loc[0, "station_nm"], "Hudson River at Albany")
        self.assertEqual(len(result), 1)

    def test_region_parameter_screen_checks_period(self):
        features = [
            {"properties": {
                "monitoring_location_id": "USGS-01358000",
                "computation_identifier": "Instantaneous",
                "begin_utc": "2010-01-01T00:00:00Z",
                "end_utc": "2020-01-01T00:00:00Z",
            }},
            {"properties": {
                "monitoring_location_id": "USGS-01359000",
                "computation_identifier": "Instantaneous",
                "begin_utc": "2000-01-01T00:00:00Z",
                "end_utc": "2010-01-01T00:00:00Z",
            }},
        ]
        bbox = (-74.0, 42.0, -73.0, 43.0)
        with patch.object(download_usgs, "iter_ogc_features", return_value=features) as request:
            matched = download_usgs.discover_station_ids_with_parameter_in_bbox(
                bbox, "00060", "2017-01-01", "2017-01-05", session=object()
            )
        self.assertEqual(matched, {"01358000"})
        self.assertEqual(request.call_args.args[2], {
            "parameter_code": "00060",
            "bbox": "-74.0,42.0,-73.0,43.0",
        })

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

    def test_rate_limit_retries_with_api_key_after_retry_after(self):
        limited = Mock(status_code=429, headers={"Retry-After": "90"})
        success = Mock(status_code=200)
        success.json.return_value = {"features": []}
        session = Mock(headers={"X-Api-Key": "test-key"})
        session.get.side_effect = [limited, success]

        with patch.object(download_usgs.time, "sleep") as sleep:
            result = download_usgs.request_json(session, "https://example.test/items")

        self.assertEqual(result, {"features": []})
        self.assertEqual(session.get.call_count, 2)
        self.assertIn("X-Api-Key", session.headers)
        sleep.assert_called_once_with(90)

    def test_rate_limit_without_retry_after_waits_for_key_limit(self):
        limited = Mock(status_code=429, headers={})
        success = Mock(status_code=200)
        success.json.return_value = {"features": []}
        session = Mock(headers={"X-Api-Key": "test-key"})
        session.get.side_effect = [limited, success]

        with patch.object(download_usgs.time, "sleep") as sleep:
            download_usgs.request_json(session, "https://example.test/items")

        self.assertIn("X-Api-Key", session.headers)
        sleep.assert_called_once_with(3600.0)

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

    def test_temperature_download_filters_unavailable_and_below_minus_two(self):
        observations = [
            {
                "geometry": {"coordinates": [-81.0, 31.0]},
                "properties": {
                    "monitoring_location_id": "USGS-02228000",
                    "time_series_id": "series-1",
                    "time": f"2017-01-01T{hour:02}:00:00Z",
                    "value": value,
                    "qualifier": qualifier,
                    "unit_of_measure": "deg C",
                },
            }
            for hour, value, qualifier in [
                (0, "-118", ["UNAVAIL"]),
                (1, "8", "UNAVAIL"),
                (2, "-2.01", []),
                (3, "-2", []),
                (4, "10", []),
            ]
        ]
        for downloader in (
            download_usgs.download_station_parameter,
            download_usgs.download_station_group_parameter,
            download_usgs.download_state_parameter,
        ):
            with self.subTest(downloader=downloader.__name__):
                with patch.object(
                    download_usgs, "iter_ogc_features", return_value=observations
                ):
                    result = downloader(
                        "02228000" if downloader is not download_usgs.download_state_parameter else "Georgia",
                        "00010", "2017-01-01", "2017-01-02", session=object()
                    )
                self.assertEqual(result["value"].tolist(), [-2.0, 10.0])

    def test_unavailable_flow_is_rejected_without_temperature_cutoff(self):
        features = [
            {
                "properties": {
                    "monitoring_location_id": "USGS-01358000",
                    "time": f"2017-01-01T{hour:02}:00:00Z",
                    "value": value,
                    "qualifier": qualifier,
                }
            }
            for hour, value, qualifier in [
                (0, "50", ["UNAVAIL"]),
                (1, "-3", []),
            ]
        ]
        with patch.object(download_usgs, "iter_ogc_features", return_value=features):
            result = download_usgs.download_station_group_parameter(
                ["01358000"], "00060", "2017-01-01", "2017-01-02",
                session=object(),
            )

        self.assertEqual(result["value"].tolist(), [-3.0])

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
