"""Characterization tests for automatic source-temperature support."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import source_temperature
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    MANUAL_NWM_TO_USGS_FLOW,
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_STATION_BY_NAME,
)


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


if __name__ == "__main__":
    unittest.main()
