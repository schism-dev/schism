"""Characterization tests for source/sink spatial helper operations."""

import unittest

import numpy as np

from stofs3d_setup.ops.Source_sink import spatial_corrections as spatial


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


if __name__ == "__main__":
    unittest.main()
