"""Characterization tests for automatic source-temperature support."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import source_temperature
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    MANUAL_NWM_TO_USGS_FLOW,
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_STATION_BY_NAME,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
)


class _Grid:
    xctr = np.array([0.0])
    yctr = np.array([0.0])

    def compute_ctr(self):
        pass


class SourceTemperatureTests(unittest.TestCase):
    def test_source_mapping_loader_normalizes_elements_and_feature_ids(self):
        with TemporaryDirectory() as directory:
            mapping_file = Path(directory) / "sources.json"
            mapping_file.write_text(
                json.dumps({"10": [1, "2"], "20": 3}),
                encoding="utf-8",
            )

            mapping = source_temperature._load_relocated_source_fids(directory)

        self.assertEqual(mapping, {10: [1, 2], 20: [3]})

    def test_source_mapping_loader_requires_sources_json(self):
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "sources.json"):
                source_temperature._load_relocated_source_fids(directory)

    def test_explicit_hudson_station_mappings_are_parameter_specific(self):
        self.assertEqual(USGS_FLOW_STATION_BY_NAME["Hudson River"], "01358000")
        self.assertEqual(
            USGS_TEMPERATURE_STATION_BY_NAME["Hudson River"],
            "01359139",
        )
        self.assertEqual(MANUAL_NWM_TO_USGS_FLOW[6186156], "01358000")

    def test_public_stage_changes_temperature_without_changing_flow(self):
        time = np.array([0.0, 3600.0])
        original = _build_source_sink(
            source_eles=[1],
            source_time_and_data=(time, np.array([[2.0], [3.0]])),
            msource_data_list=[
                (time, np.array([[10.0], [11.0]])),
                (time, np.array([[0.0], [0.0]])),
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )

        def replace_temperature(**kwargs):
            tracers = kwargs["msource_data_list"]
            tracers[0][1][:, 0] = [20.0, 21.0]
            return tracers, 1

        with patch.object(
            source_temperature,
            "_replace_all_relocated_source_temperatures",
            side_effect=replace_temperature,
        ):
            corrected, count = (
                source_temperature.replace_source_temperatures_with_usgs(
                    base_ss=original,
                    hgrid=_Grid(),
                    source_mapping_dir="unused",
                    start_time="2020-01-01",
                    usgs_cache_folder="unused",
                    nwm_shapefile="unused",
                )
            )

        self.assertEqual(count, 1)
        self.assertEqual(np.asarray(corrected.source_eles).tolist(), [1])
        np.testing.assert_array_equal(
            corrected.vsource.df.values,
            [[2.0], [3.0]],
        )
        np.testing.assert_array_equal(
            corrected.msource[0].df.values,
            [[20.0], [21.0]],
        )
        np.testing.assert_array_equal(
            original.msource[0].df.values,
            [[10.0], [11.0]],
        )


if __name__ == "__main__":
    unittest.main()
