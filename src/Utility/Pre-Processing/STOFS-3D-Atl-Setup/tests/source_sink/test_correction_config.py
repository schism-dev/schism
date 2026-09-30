"""Characterization tests for source/sink correction configuration."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stofs3d_setup.ops.Source_sink import correction_config


class CorrectionConfigTests(unittest.TestCase):
    def test_replace_override_defaults_are_normalized(self):
        points = correction_config._normalize_replace_relocated_points(
            {
                "replace_only_source_locations": [
                    {
                        "name": "Delaware",
                        "x": -74.9,
                        "y": 40.3,
                    }
                ]
            }
        )

        self.assertEqual(points[0]["max_search_radius_m"], 500.0)
        self.assertTrue(points[0]["replace_flow"])
        self.assertTrue(points[0]["replace_temperature"])
        self.assertFalse(points[0]["allow_negative_sink"])

    def test_exclusion_remove_string_becomes_list(self):
        points = correction_config._normalize_exclude_points(
            {
                "exclude_source_sink_locations": [
                    {
                        "name": "test",
                        "x": 1,
                        "y": 2,
                        "remove": "source",
                    }
                ]
            }
        )

        self.assertEqual(points[0]["remove"], ["source"])
        self.assertEqual(points[0]["radius_m"], 500.0)

    def test_relative_region_path_resolves_from_yaml_directory(self):
        with TemporaryDirectory() as directory:
            directory = Path(directory)
            region_file = directory / "region.rgn"
            region_file.touch()

            regions = correction_config._normalize_zero_source_regions(
                {"zero_source_regions": ["region.rgn"]},
                yaml_dir=directory,
            )

        self.assertEqual(regions[0]["region_file"], region_file.resolve())


if __name__ == "__main__":
    unittest.main()
