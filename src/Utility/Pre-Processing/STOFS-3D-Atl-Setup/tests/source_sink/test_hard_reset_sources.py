"""Characterization tests for explicit USGS source overrides."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import hard_reset_sources
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
)


class _IdentityTransformer:
    def transform(self, x, y):
        return x, y


class _Grid:
    xctr = np.array([0.0])
    yctr = np.array([0.0])

    def compute_ctr(self):
        pass


class _RouteGrid:
    xctr = np.array([-74.0])
    yctr = np.array([40.0])

    def compute_ctr(self):
        pass


class SourceOverrideTests(unittest.TestCase):
    def test_optional_wrapper_applies_flow_then_temperature(self):
        original = object()
        flow_corrected = object()
        temperature_corrected = object()
        points = [{"name": "Delaware"}]
        calls = []

        def apply_flow(**kwargs):
            calls.append("flow")
            self.assertIs(kwargs["base_ss"], original)
            self.assertEqual(kwargs["points"], points)
            return flow_corrected, 1

        def apply_temperature(**kwargs):
            calls.append("temperature")
            self.assertIs(kwargs["base_ss"], flow_corrected)
            self.assertEqual(kwargs["points"], points)
            return temperature_corrected, 1

        with patch.object(
            hard_reset_sources,
            "load_source_override_points",
            return_value=points,
        ), patch.object(
            hard_reset_sources,
            "apply_source_flow_overrides",
            side_effect=apply_flow,
        ), patch.object(
            hard_reset_sources,
            "apply_source_temperature_overrides",
            side_effect=apply_temperature,
        ):
            result = hard_reset_sources.apply_hard_reset_sources_with_usgs(
                base_ss=original,
                hgrid=_Grid(),
                enabled=True,
                correction_info="overrides.yml",
                start_time="2020-01-01",
                usgs_cache_folder="cache",
            )

        self.assertEqual(calls, ["flow", "temperature"])
        self.assertEqual(result, (temperature_corrected, 1, 1))

    def test_optional_wrapper_skips_disabled_overrides(self):
        original = object()

        with patch.object(
            hard_reset_sources, "load_source_override_points"
        ) as load_points:
            result = hard_reset_sources.apply_hard_reset_sources_with_usgs(
                base_ss=original,
                hgrid=_Grid(),
                enabled=False,
                correction_info=None,
                start_time="2020-01-01",
                usgs_cache_folder=None,
            )

        load_points.assert_not_called()
        self.assertEqual(result, (original, 0, 0))

    def test_flow_override_replaces_supported_records_directly(self):
        model_time = np.array([0.0, 3600.0])
        observations = pd.Series(
            [10.0, 20.0],
            index=pd.to_datetime(
                ["2020-01-01T00:00Z", "2020-01-01T01:00Z"]
            ),
        )
        point = {
            "name": "Delaware",
            "x": 0.0,
            "y": 0.0,
            "flow_station_id": "01463500",
            "temperature_station_id": "01463500",
            "max_search_radius_m": 500.0,
            "replace_flow": True,
            "replace_temperature": False,
            "allow_negative_sink": False,
        }

        with patch.object(
            hard_reset_sources,
            "_get_usgs_flow",
            return_value=(observations, "01463500"),
        ) as get_flow:
            result = hard_reset_sources._replace_existing_relocated_source(
                point=point,
                source_eles=[10],
                source_time_and_data=(
                    model_time,
                    np.array([[1.0], [2.0]]),
                ),
                msource_data_list=[],
                sink_eles=[],
                sink_time_and_data=None,
                xctr=np.zeros(10),
                yctr=np.zeros(10),
                transformer=_IdentityTransformer(),
                start_time=pd.Timestamp("2020-01-01T00:00Z"),
                usgs_cache_folder=None,
            )

        self.assertNotIn("download_policy", get_flow.call_args.kwargs)

        np.testing.assert_allclose(
            result[0][1][:, 0],
            np.array([10.0, 20.0]) * CFS_TO_CMS,
        )

    def test_override_requires_source_within_configured_radius(self):
        point = {
            "name": "Delaware",
            "x": 1000.0,
            "y": 0.0,
            "flow_station_id": "01463500",
            "temperature_station_id": "01463500",
            "max_search_radius_m": 10.0,
            "replace_flow": True,
            "replace_temperature": False,
            "allow_negative_sink": False,
        }

        with self.assertRaisesRegex(ValueError, "no relocated source found"):
            hard_reset_sources._replace_existing_relocated_source(
                point=point,
                source_eles=[1],
                source_time_and_data=(
                    np.array([0.0]),
                    np.array([[1.0]]),
                ),
                msource_data_list=[],
                sink_eles=[],
                sink_time_and_data=None,
                xctr=np.array([0.0]),
                yctr=np.array([0.0]),
                transformer=_IdentityTransformer(),
                start_time=pd.Timestamp("2020-01-01T00:00Z"),
                usgs_cache_folder=None,
            )

    def test_hard_reset_follows_station_route_to_first_grid_source(self):
        network = gpd.GeoDataFrame(
            {"to": [6186164, 6189262, 6189262]},
            geometry=[
                LineString([(-74.02, 40.0), (-74.01, 40.0)]),
                LineString([(-74.01, 40.0), (-74.0, 40.0)]),
                LineString([(-74.0, 40.0), (-73.99, 40.0)]),
            ],
            index=[6186156, 6186164, 6189262],
            crs="EPSG:4326",
        )
        base_ss = type("SourceSink", (), {"source_eles": [1]})()

        with TemporaryDirectory() as directory:
            source_map = Path(directory) / "sources.json"
            source_map.write_text('{"1": [6189262]}', encoding="utf-8")
            with patch.object(
                hard_reset_sources,
                "preprocess_nwm_shp",
                return_value=network,
            ):
                resolved = (
                    hard_reset_sources.resolve_hard_reset_source_targets(
                        base_ss=base_ss,
                        hgrid=_RouteGrid(),
                        points=[
                            {
                                "name": "Hudson River",
                                "x": -73.7,
                                "y": 42.7,
                                "flow_station_id": "01358000",
                                "temperature_station_id": "01359139",
                                "max_search_radius_m": 2500.0,
                                "replace_flow": True,
                                "replace_temperature": True,
                            }
                        ],
                        source_mapping_dir=directory,
                    )
                )

        self.assertEqual(resolved[0]["_target_element_id"], 1)
        self.assertEqual(resolved[0]["_target_nwm_feature_id"], 6189262)
        self.assertLessEqual(resolved[0]["_target_distance_m"], 2500.0)

    def test_unmapped_station_is_associated_from_usgs_location(self):
        station_id = "00000001"
        network = gpd.GeoDataFrame(
            {"to": [9001]},
            geometry=[LineString([(-74.01, 40.0), (-73.99, 40.0)])],
            index=[9001],
            crs="EPSG:4326",
        )
        locations = pd.DataFrame(
            [{"station_id": station_id, "lon": -74.0, "lat": 40.0}]
        )

        with (
            patch.object(
                hard_reset_sources,
                "NWM_TO_USGS_FLOW_ADJUSTMENT",
                {},
            ),
            patch.object(
                hard_reset_sources,
                "NWM_TO_USGS_TEMPERATURE_SEARCH",
                {},
            ),
            patch.object(
                hard_reset_sources,
                "discover_monitoring_locations",
                return_value=locations,
            ),
        ):
            result = hard_reset_sources._station_feature_id(
                {
                    "replace_flow": True,
                    "flow_station_id": station_id,
                },
                network,
            )

        self.assertEqual(result, (station_id, 9001))

    def test_downstream_source_match_keeps_radius_against_flowline(self):
        network = gpd.GeoDataFrame(
            {"to": [6189262, 6189262]},
            geometry=[
                LineString([(-74.02, 40.0), (-74.01, 40.0)]),
                LineString([(-74.01, 40.0), (-74.0, 40.0)]),
            ],
            index=[6186156, 6189262],
            crs="EPSG:4326",
        )
        xctr = np.array([-74.0])
        yctr = np.array([40.1])

        with self.assertRaisesRegex(
            ValueError,
            "first downstream grid source .* outside the 2500.0 m search radius",
        ):
            hard_reset_sources._search_downstream_source(
                start_feature_id=6186156,
                nwm_network=network,
                source_feature_mapping={6189262: [1]},
                source_element_ids={1},
                xctr=xctr,
                yctr=yctr,
                transformer=hard_reset_sources._make_transformer(),
                max_search_radius_m=2500.0,
            )

    def test_public_stage_processes_points_in_order(self):
        time = np.array([0.0])
        original = _build_source_sink(
            source_eles=[1],
            source_time_and_data=(time, np.array([[1.0]])),
            msource_data_list=[
                (time, np.array([[10.0]])),
                (time, np.array([[0.0]])),
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )
        processed = []

        def record_override(**kwargs):
            processed.append(kwargs["point"]["name"])
            source_time, source_data = kwargs["source_time_and_data"]
            source_data[:, 0] += 1.0
            return (
                (source_time, source_data),
                kwargs["msource_data_list"],
                kwargs["sink_eles"],
                kwargs["sink_time_and_data"],
            )

        with patch.object(
            hard_reset_sources,
            "_replace_existing_relocated_source",
            side_effect=record_override,
        ):
            corrected, count = hard_reset_sources.apply_source_overrides(
                base_ss=original,
                hgrid=_Grid(),
                points=[{"name": "first"}, {"name": "second"}],
                start_time="2020-01-01",
                usgs_cache_folder=None,
            )

        self.assertEqual(count, 2)
        self.assertEqual(processed, ["first", "second"])
        np.testing.assert_array_equal(corrected.vsource.df.values, [[3.0]])
        np.testing.assert_array_equal(original.vsource.df.values, [[1.0]])

    def test_flow_stage_does_not_download_or_change_temperature(self):
        time = np.array([0.0, 3600.0])
        original = _build_source_sink(
            source_eles=[1],
            source_time_and_data=(time, np.array([[1.0], [2.0]])),
            msource_data_list=[
                (time, np.array([[10.0], [11.0]])),
                (time, np.array([[0.0], [0.0]])),
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )
        point = {
            "name": "Delaware",
            "x": 0.0,
            "y": 0.0,
            "flow_station_id": "01463500",
            "temperature_station_id": "01463500",
            "max_search_radius_m": 500.0,
            "replace_flow": True,
            "replace_temperature": True,
            "allow_negative_sink": False,
        }
        observations = pd.Series(
            [10.0, 20.0],
            index=pd.to_datetime(
                ["2020-01-01T00:00Z", "2020-01-01T01:00Z"]
            ),
        )

        with (
            patch.object(
                hard_reset_sources,
                "_get_usgs_flow",
                return_value=(observations, "01463500"),
            ),
            patch.object(hard_reset_sources, "_get_usgs_temperature") as get_temp,
        ):
            corrected, count = hard_reset_sources.apply_source_flow_overrides(
                base_ss=original,
                hgrid=_Grid(),
                points=[point],
                start_time=pd.Timestamp("2020-01-01T00:00Z"),
                usgs_cache_folder=None,
            )

        self.assertEqual(count, 1)
        get_temp.assert_not_called()
        np.testing.assert_allclose(
            corrected.vsource.df.values[:, 0],
            np.array([10.0, 20.0]) * CFS_TO_CMS,
        )
        np.testing.assert_array_equal(
            corrected.msource[0].df.values,
            original.msource[0].df.values,
        )
        np.testing.assert_array_equal(
            corrected.msource[1].df.values,
            original.msource[1].df.values,
        )

    def test_temperature_override_rejects_partial_series_in_both_modes(self):
        time = np.array([0.0, 3600.0, 7200.0])
        original = _build_source_sink(
            source_eles=[1],
            source_time_and_data=(time, np.ones((3, 1))),
            msource_data_list=[
                (time, np.full((3, 1), -9999.0)),
                (time, np.zeros((3, 1))),
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )
        point = {
            "name": "Delaware",
            "x": 0.0,
            "y": 0.0,
            "flow_station_id": "01463500",
            "temperature_station_id": "01463500",
            "max_search_radius_m": 500.0,
            "replace_flow": False,
            "replace_temperature": True,
            "allow_negative_sink": False,
        }
        partial_temperature = pd.Series(
            [12.0],
            index=pd.to_datetime(["2020-01-01T00:00Z"]),
        )

        for mode in ("first_usable", "discharge_weighted"):
            with self.subTest(mode=mode), patch.object(
                hard_reset_sources,
                "_get_usgs_temperature",
                return_value=(partial_temperature, "01463500"),
            ):
                corrected, count = (
                    hard_reset_sources.apply_source_temperature_overrides(
                        base_ss=original,
                        hgrid=_Grid(),
                        points=[point],
                        start_time=pd.Timestamp("2020-01-01T00:00Z"),
                        usgs_cache_folder=None,
                        temperature_pooling=mode,
                    )
                )
                self.assertEqual(count, 1)
                np.testing.assert_array_equal(
                    corrected.msource[0].df.values,
                    original.msource[0].df.values,
                )
                np.testing.assert_array_equal(
                    corrected.vsource.df.values,
                    original.vsource.df.values,
                )


if __name__ == "__main__":
    unittest.main()
