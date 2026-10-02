"""Regression tests for combining source/sink objects."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from pylib_experimental.schism_file import TimeHistory, source_sink


class SourceSinkAdditionTests(unittest.TestCase):
    def test_addition_is_independent_and_roundtrips_through_netcdf(self):
        time = np.array([0.0, 3600.0])
        base = source_sink(
            vsource=TimeHistory(
                data_array=np.c_[time, [[1.0], [2.0]]],
                columns=[2],
            ),
            vsink=TimeHistory(
                data_array=np.c_[time, [[-0.5], [-1.0]]],
                columns=[4],
            ),
            msource=[
                TimeHistory(
                    data_array=np.c_[time, [[10.0], [20.0]]],
                    columns=[2],
                )
            ],
        )
        background = source_sink(
            vsource=TimeHistory(
                data_array=np.c_[time, [[3.0], [4.0]]],
                columns=[3],
            ),
            vsink=TimeHistory(
                data_array=np.c_[time, [[-1.5], [-2.0]]],
                columns=[5],
            ),
            msource=[
                TimeHistory(
                    data_array=np.c_[time, [[30.0], [40.0]]],
                    columns=[3],
                )
            ],
        )

        combined = base + background

        np.testing.assert_array_equal(combined.source_eles, [2, 3])
        np.testing.assert_array_equal(combined.sink_eles, [4, 5])
        np.testing.assert_array_equal(
            combined.vsource.data,
            [[1.0, 3.0], [2.0, 4.0]],
        )
        np.testing.assert_array_equal(
            combined.msource[0].data,
            [[10.0, 30.0], [20.0, 40.0]],
        )
        np.testing.assert_array_equal(
            combined.vsink.data,
            [[-0.5, -1.5], [-1.0, -2.0]],
        )

        with tempfile.TemporaryDirectory() as output_dir:
            combined.writer(output_dir)
            restored = source_sink.from_ncfile(Path(output_dir) / "source.nc")

        np.testing.assert_array_equal(restored.source_eles, [2, 3])
        np.testing.assert_array_equal(restored.sink_eles, [4, 5])
        np.testing.assert_allclose(restored.vsource.data, combined.vsource.data)
        np.testing.assert_allclose(
            restored.msource[0].data, combined.msource[0].data
        )
        np.testing.assert_allclose(restored.vsink.data, combined.vsink.data)

        combined.vsource.df.iloc[0, 0] = 99.0
        combined.msource[0].df.iloc[0, 0] = 99.0
        combined.vsink.df.iloc[0, 0] = -99.0
        np.testing.assert_array_equal(base.vsource.data[:, 0], [1.0, 2.0])
        np.testing.assert_array_equal(
            base.msource[0].data[:, 0], [10.0, 20.0]
        )
        np.testing.assert_array_equal(base.vsink.data[:, 0], [-0.5, -1.0])
        np.testing.assert_array_equal(background.vsource.data[:, 0], [3.0, 4.0])
        np.testing.assert_array_equal(
            background.msource[0].data[:, 0], [30.0, 40.0]
        )
        np.testing.assert_array_equal(
            background.vsink.data[:, 0], [-1.5, -2.0]
        )


if __name__ == "__main__":
    unittest.main()
