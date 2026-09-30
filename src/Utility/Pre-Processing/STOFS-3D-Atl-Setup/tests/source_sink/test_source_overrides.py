"""Characterization tests for explicit USGS source overrides."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import source_overrides
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
)


class _IdentityTransformer:
    def transform(self, x, y):
        return x, y


class SourceOverrideTests(unittest.TestCase):
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
        }

        with patch.object(
            source_overrides,
            "_get_usgs_forcing",
            return_value=(observations, None, "01463500", "01463500"),
        ):
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


if __name__ == "__main__":
    unittest.main()
