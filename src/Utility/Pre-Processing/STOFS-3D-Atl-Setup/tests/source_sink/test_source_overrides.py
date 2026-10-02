"""Characterization tests for explicit USGS source overrides."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import source_overrides
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
            source_overrides,
            "load_source_override_points",
            return_value=points,
        ), patch.object(
            source_overrides,
            "apply_source_flow_overrides",
            side_effect=apply_flow,
        ), patch.object(
            source_overrides,
            "apply_source_temperature_overrides",
            side_effect=apply_temperature,
        ):
            result = source_overrides.apply_selected_source_overrides(
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
            source_overrides, "load_source_override_points"
        ) as load_points:
            result = source_overrides.apply_selected_source_overrides(
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
            "max_search_radius_m": 500.0,
            "replace_flow": True,
            "replace_temperature": False,
            "allow_negative_sink": False,
            "usgs_download": {
                "primary_chunk_days": 100,
                "retry_chunk_days": [20, 10, 5],
            },
        }

        with patch.object(
            source_overrides,
            "_get_usgs_flow",
            return_value=(observations, "01463500"),
        ) as get_flow:
            result = source_overrides._replace_existing_relocated_source(
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

        self.assertEqual(
            get_flow.call_args.kwargs["download_policy"],
            point["usgs_download"],
        )

        np.testing.assert_allclose(
            result[0][1][:, 0],
            np.array([10.0, 20.0]) * CFS_TO_CMS,
        )

    def test_override_requires_source_within_configured_radius(self):
        point = {
            "name": "Delaware",
            "x": 1000.0,
            "y": 0.0,
            "max_search_radius_m": 10.0,
            "replace_flow": True,
            "replace_temperature": False,
            "allow_negative_sink": False,
        }

        with self.assertRaisesRegex(ValueError, "no relocated source found"):
            source_overrides._replace_existing_relocated_source(
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
            source_overrides,
            "_replace_existing_relocated_source",
            side_effect=record_override,
        ):
            corrected, count = source_overrides.apply_source_overrides(
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
                source_overrides,
                "_get_usgs_flow",
                return_value=(observations, "01463500"),
            ),
            patch.object(source_overrides, "_get_usgs_temperature") as get_temp,
        ):
            corrected, count = source_overrides.apply_source_flow_overrides(
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

    def test_weighted_mode_rejects_partial_selected_temperature_override(self):
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
            "max_search_radius_m": 500.0,
            "replace_flow": False,
            "replace_temperature": True,
            "allow_negative_sink": False,
        }
        partial_temperature = pd.Series(
            [12.0],
            index=pd.to_datetime(["2020-01-01T00:00Z"]),
        )

        with patch.object(
            source_overrides,
            "_get_usgs_temperature",
            return_value=(partial_temperature, "01463500"),
        ):
            corrected, count = (
                source_overrides.apply_source_temperature_overrides(
                    base_ss=original,
                    hgrid=_Grid(),
                    points=[point],
                    start_time=pd.Timestamp("2020-01-01T00:00Z"),
                    usgs_cache_folder=None,
                    temperature_pooling="discharge_weighted",
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
