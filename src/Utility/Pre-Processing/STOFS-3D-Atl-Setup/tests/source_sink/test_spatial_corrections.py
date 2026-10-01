"""Characterization tests for source/sink spatial helper operations."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from stofs3d_setup.ops.Source_sink import spatial_corrections as spatial
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
)


class _IdentityTransformer:
    def transform(self, x, y):
        return x, y


class _Grid:
    def __init__(self):
        self.xctr = np.array([-1.0])
        self.yctr = np.array([-2.0])

    def compute_ctr(self):
        self.xctr = np.array([1.0, 2.0])
        self.yctr = np.array([3.0, 4.0])


class SpatialCorrectionTests(unittest.TestCase):
    def test_configured_region_stage_owns_loading_and_path_resolution(self):
        with TemporaryDirectory() as directory:
            directory = Path(directory)
            region_file = directory / "test.rgn"
            region_file.touch()
            config_file = directory / "regions.yml"
            config_file.write_text(
                "zero_source_regions:\n"
                "  - name: test region\n"
                "    region_file: test.rgn\n",
                encoding="utf-8",
            )

            with patch.object(
                spatial,
                "zero_sources_in_regions",
                return_value=("corrected", 1),
            ) as zero_regions:
                result = spatial.zero_configured_source_regions(
                    base_ss="base",
                    hgrid="grid",
                    correction_info=config_file,
                )

        self.assertEqual(result, ("corrected", 1))
        self.assertEqual(
            zero_regions.call_args.kwargs["regions"],
            [{"name": "test region", "region_file": region_file}],
        )

    def test_compute_grid_centers_restores_existing_grid_attributes(self):
        grid = _Grid()

        xctr, yctr = spatial._compute_grid_centers(grid)

        np.testing.assert_array_equal(xctr, [1.0, 2.0])
        np.testing.assert_array_equal(yctr, [3.0, 4.0])
        np.testing.assert_array_equal(grid.xctr, [-1.0])
        np.testing.assert_array_equal(grid.yctr, [-2.0])

    def test_nearest_element_searches_only_candidate_elements(self):
        element, distance = spatial._nearest_element(
            x=3.4,
            y=0.0,
            candidate_element_ids=[1, 2],
            xctr=np.array([0.0, 3.0, 3.3]),
            yctr=np.zeros(3),
            transformer=_IdentityTransformer(),
        )

        self.assertEqual(element, 2)
        self.assertAlmostEqual(distance, 0.4)

    def test_elements_within_radius_returns_only_candidates_in_radius(self):
        elements = spatial._elements_within_radius(
            x=2.0,
            y=0.0,
            radius_m=2.1,
            candidate_element_ids=[1, 2, 4],
            xctr=np.array([0.0, 3.0, 2.0, 10.0]),
            yctr=np.zeros(4),
            transformer=_IdentityTransformer(),
        )

        self.assertEqual(set(elements), {1, 2})

    def test_empty_candidate_list_has_no_match(self):
        element, distance = spatial._nearest_element(
            x=0.0,
            y=0.0,
            candidate_element_ids=[],
            xctr=np.array([]),
            yctr=np.array([]),
            transformer=_IdentityTransformer(),
        )

        self.assertIsNone(element)
        self.assertTrue(np.isinf(distance))

    def test_public_region_stage_preserves_tracers_and_element_order(self):
        time = np.array([0.0, 3600.0])
        original = _build_source_sink(
            source_eles=[1, 2],
            source_time_and_data=(
                time,
                np.array([[1.0, 2.0], [3.0, 4.0]]),
            ),
            msource_data_list=[
                (time, np.array([[10.0, 20.0], [30.0, 40.0]]))
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )

        def zero_first_source(**kwargs):
            source_time, source_data = kwargs["source_time_and_data"]
            source_data[:, 0] = 0.0
            return (source_time, source_data), 1

        with patch.object(
            spatial,
            "_zero_sources_inside_regions",
            side_effect=zero_first_source,
        ):
            corrected, count = spatial.zero_sources_in_regions(
                original,
                _Grid(),
                regions=[{"name": "test"}],
            )

        self.assertEqual(count, 1)
        self.assertEqual(np.asarray(corrected.source_eles).tolist(), [1, 2])
        np.testing.assert_array_equal(
            corrected.vsource.df.values,
            [[0.0, 2.0], [0.0, 4.0]],
        )
        np.testing.assert_array_equal(
            corrected.msource[0].df.values,
            [[10.0, 20.0], [30.0, 40.0]],
        )
        np.testing.assert_array_equal(
            original.vsource.df.values,
            [[1.0, 2.0], [3.0, 4.0]],
        )


if __name__ == "__main__":
    unittest.main()
