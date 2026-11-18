"""Tests for background-sink post-processing."""

import unittest

import numpy as np

from pylib_experimental.schism_file import source_sink, TimeHistory
from stofs3d_setup.ops.Source_sink.Constant_sinks.background_sink import (
    remove_overlapping_background_sinks,
)


class BackgroundSinkTests(unittest.TestCase):
    def test_sinks_overlapping_base_sources_are_removed(self):
        time = np.array([0.0, 3600.0])
        base = source_sink(
            vsource=TimeHistory(
                data_array=np.c_[time, [[1.0], [1.0]]],
                columns=["2"],
            ),
            vsink=None,
            msource=[
                TimeHistory(
                    data_array=np.c_[time, [[-9999.0], [-9999.0]]],
                    columns=["2"],
                )
            ],
        )
        background = source_sink(
            vsource=None,
            vsink=TimeHistory(
                data_array=np.c_[time, [[-1.0, -2.0], [-1.0, -2.0]]],
                columns=["2", "3"],
            ),
            msource=None,
        )

        result = remove_overlapping_background_sinks(background, base)

        np.testing.assert_array_equal(result.sink_eles, [3])
        np.testing.assert_array_equal(result.vsink.data[:, 0], [-2.0, -2.0])


if __name__ == "__main__":
    unittest.main()
