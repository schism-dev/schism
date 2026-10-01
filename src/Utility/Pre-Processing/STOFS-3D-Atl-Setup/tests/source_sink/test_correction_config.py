"""Characterization tests for source/sink correction configuration."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stofs3d_setup.ops.Source_sink import correction_config
from stofs3d_setup.ops.Source_sink.Replace_with_USGS import override_config
from stofs3d_setup.ops.Source_sink.Spatial_corrections import region_config
from stofs3d_setup.config.stofs3d_atl_config import ConfigStofs3dAtlantic


class CorrectionConfigTests(unittest.TestCase):
    def test_independent_correction_stages_are_disabled_by_default(self):
        config = ConfigStofs3dAtlantic()

        self.assertFalse(config.replace_source_temperature_with_usgs)
        self.assertEqual(config.source_temperature_pooling, "first_usable")
        self.assertFalse(config.replace_selected_sources_with_usgs)
        self.assertFalse(config.zero_configured_source_regions)
        self.assertIsNone(config.selected_source_override_info)
        self.assertIsNone(config.zero_source_region_info)

    def test_v7p4_enables_preferred_correction_stages(self):
        config = ConfigStofs3dAtlantic.v7p4()

        self.assertTrue(config.replace_source_temperature_with_usgs)
        self.assertEqual(config.source_temperature_pooling, "first_usable")
        self.assertFalse(config.replace_selected_sources_with_usgs)
        self.assertFalse(config.zero_configured_source_regions)
        self.assertNotEqual(
            config.selected_source_override_info,
            config.artificial_island_source_sink_info,
        )
        self.assertNotEqual(
            config.zero_source_region_info,
            config.artificial_island_source_sink_info,
        )
        self.assertNotEqual(
            config.selected_source_override_info,
            config.zero_source_region_info,
        )

    def test_v7p4_stage_files_contain_only_owned_sections(self):
        config = ConfigStofs3dAtlantic.v7p4()

        overrides = correction_config.load_source_sink_corrections(
            config.selected_source_override_info
        )
        regions = correction_config.load_source_sink_corrections(
            config.zero_source_region_info
        )
        islands = correction_config.load_source_sink_corrections(
            config.artificial_island_source_sink_info
        )

        self.assertEqual(set(overrides), {"replace_only_source_locations"})
        self.assertEqual(set(regions), {"zero_source_regions"})
        self.assertEqual(
            set(islands),
            {
                "remove_source_locations_in_artificial_island",
                "force_source_sink_locations",
                "large_constant_sink_artificial_island_locations",
            },
        )

        override_points = override_config.source_override_points(overrides)
        self.assertEqual(
            [point["name"] for point in override_points],
            ["Delaware", "Hudson River"],
        )
        self.assertEqual(
            override_points[1]["usgs_download"]["retry_chunk_days"],
            [20, 10, 5],
        )
        region_points = region_config.zero_source_regions(
            regions,
            config_dir=config.zero_source_region_info.parent,
        )
        self.assertEqual(region_points[0]["name"], "Upstream Savannah")
        self.assertTrue(region_points[0]["region_file"].is_file())

    def test_replace_override_defaults_are_normalized(self):
        points = override_config._normalize_replace_relocated_points(
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

    def test_relative_region_path_resolves_from_yaml_directory(self):
        with TemporaryDirectory() as directory:
            directory = Path(directory)
            region_file = directory / "region.rgn"
            region_file.touch()

            regions = region_config._normalize_zero_source_regions(
                {"zero_source_regions": ["region.rgn"]},
                yaml_dir=directory,
            )

        self.assertEqual(regions[0]["region_file"], region_file.resolve())

if __name__ == "__main__":
    unittest.main()
