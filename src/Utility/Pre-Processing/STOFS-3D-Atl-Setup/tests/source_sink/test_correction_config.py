"""Characterization tests for source/sink correction configuration."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stofs3d_setup.ops.Source_sink import correction_config
from stofs3d_setup.ops.Source_sink.Replace_with_USGS import override_config
from stofs3d_setup.ops.Source_sink.Patch_artificial_island import island_config
from stofs3d_setup.ops.Source_sink.Spatial_corrections import region_config
from stofs3d_setup.config.stofs3d_atl_config import ConfigStofs3dAtlantic


class CorrectionConfigTests(unittest.TestCase):
    def test_independent_correction_stages_are_disabled_by_default(self):
        config = ConfigStofs3dAtlantic()

        self.assertFalse(config.replace_source_temperature_with_usgs)
        self.assertEqual(config.source_temperature_pooling, "first_usable")
        self.assertFalse(config.hard_reset_sources_with_usgs)
        self.assertFalse(config.zero_configured_source_regions)
        self.assertIsNone(config.hard_reset_sources_with_usgs_info)
        self.assertIsNone(config.zero_source_region_info)

    def test_v7p4_enables_preferred_correction_stages(self):
        config = ConfigStofs3dAtlantic.v7p4()

        self.assertTrue(config.replace_source_temperature_with_usgs)
        self.assertEqual(config.source_temperature_pooling, "first_usable")
        self.assertFalse(config.hard_reset_sources_with_usgs)
        self.assertTrue(config.zero_configured_source_regions)
        self.assertIsNone(config.hard_reset_sources_with_usgs_info)
        self.assertIsNone(config.zero_source_region_info)

    def test_v7p4_stage_files_contain_only_owned_sections(self):
        config = ConfigStofs3dAtlantic.v7p4()
        override_info = (
            Path(override_config.__file__).with_name(
                "hard_reset_sources_with_usgs.yml"
            )
        )
        region_info = Path(region_config.__file__).with_name(
            "zero_source_regions.yml"
        )

        overrides = correction_config.load_source_sink_corrections(
            override_info
        )
        regions = correction_config.load_source_sink_corrections(
            region_info
        )
        islands = correction_config.load_source_sink_corrections(
            config.artificial_island_source_sink_info
        )

        self.assertEqual(set(overrides), {"hard_reset_source_locations"})
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
        self.assertEqual(override_points[0]["flow_station_id"], "01463500")
        self.assertEqual(
            override_points[0]["temperature_station_id"], "01463500"
        )
        self.assertEqual(override_points[1]["flow_station_id"], "01358000")
        self.assertEqual(
            override_points[1]["temperature_station_id"], "01359139"
        )
        self.assertEqual(
            override_points[1]["usgs_download"]["retry_chunk_days"],
            [20, 10, 5],
        )
        island_points = island_config._normalize_force_points(islands)
        island_sources = {
            point["name"]: point
            for point in island_points
            if point["use_usgs_obs"]
        }
        self.assertEqual(
            island_sources["Turkey"]["flow_station_id"], "02172035"
        )
        self.assertIsNone(
            island_sources["Turkey"]["temperature_station_id"]
        )
        self.assertEqual(
            island_sources["Buffalo Bluff"]["temperature_station_id"],
            "02244040",
        )
        self.assertEqual(
            island_sources["Dunns Creek"]["flow_station_id"], "02244440"
        )
        region_points = region_config.zero_source_regions(
            regions,
            config_dir=region_info.parent,
        )
        self.assertEqual(region_points[0]["name"], "Upstream Savannah")
        self.assertTrue(region_points[0]["region_file"].is_file())

    def test_replace_override_defaults_are_normalized(self):
        points = override_config._normalize_replace_relocated_points(
            {
                "hard_reset_source_locations": [
                    {
                        "name": "Delaware",
                        "x": -74.9,
                        "y": 40.3,
                        "flow_station_id": "01463500",
                        "temperature_station_id": "01463500",
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
