"""Characterization tests for aligned source/sink component operations."""

import unittest

import numpy as np

from stofs3d_setup.ops.Source_sink import source_sink_components as components


class SourceSinkComponentTests(unittest.TestCase):
    def test_remove_source_columns_keeps_all_tracers_aligned(self):
        source_eles = [10, 20, 30]
        time = np.array([0.0, 3600.0])
        source_data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        temperature = np.array([[11.0, 12.0, 13.0], [14.0, 15.0, 16.0]])
        salinity = np.array([[21.0, 22.0, 23.0], [24.0, 25.0, 26.0]])

        result = components._remove_source_columns(
            source_eles=source_eles,
            source_time_and_data=(time, source_data),
            msource_data_list=[
                (time, temperature),
                (time, salinity),
            ],
            remove_eles={20},
        )

        kept_eles, (_, kept_source), kept_tracers = result
        self.assertEqual(kept_eles, [10, 30])
        np.testing.assert_array_equal(kept_source, source_data[:, [0, 2]])
        np.testing.assert_array_equal(
            kept_tracers[0][1], temperature[:, [0, 2]]
        )
        np.testing.assert_array_equal(
            kept_tracers[1][1], salinity[:, [0, 2]]
        )

    def test_append_source_column_keeps_tracers_aligned(self):
        time = np.array([0.0, 3600.0])
        source_eles, source_values, tracers = components._append_source_column(
            source_eles=[10],
            source_time_and_data=(time, np.array([[1.0], [2.0]])),
            msource_data_list=[(time, np.array([[11.0], [12.0]]))],
            target_ele=20,
            source_time=time,
            flow=np.array([3.0, 4.0]),
            tracer_columns=[np.array([13.0, 14.0])],
        )

        self.assertEqual(source_eles, [10, 20])
        np.testing.assert_array_equal(
            source_values[1],
            [[1.0, 3.0], [2.0, 4.0]],
        )
        np.testing.assert_array_equal(
            tracers[0][1],
            [[11.0, 13.0], [12.0, 14.0]],
        )

    def test_append_source_rejects_duplicate_element(self):
        time = np.array([0.0])
        with self.assertRaisesRegex(ValueError, "already contains a source"):
            components._append_source_column(
                source_eles=[10],
                source_time_and_data=(time, np.array([[1.0]])),
                msource_data_list=[],
                target_ele=10,
                source_time=time,
                flow=np.array([2.0]),
                tracer_columns=[],
            )

    def test_negative_source_values_move_to_sink(self):
        source_data = np.array([[-2.0], [3.0], [-1.0]])
        time = np.array([0.0, 3600.0, 7200.0])

        sink_eles, sink_values, count, minimum = (
            components._add_negative_source_part_to_sink(
                target_ele=20,
                source_column_idx=0,
                source_data=source_data,
                source_time=time,
                sink_eles=[],
                sink_time_and_data=None,
            )
        )

        self.assertEqual(sink_eles, [20])
        self.assertEqual(count, 2)
        self.assertEqual(minimum, -2.0)
        np.testing.assert_array_equal(source_data[:, 0], [0.0, 3.0, 0.0])
        np.testing.assert_array_equal(sink_values[1][:, 0], [-2.0, 0.0, -1.0])


if __name__ == "__main__":
    unittest.main()
