"""Characterization tests for USGS series handling used by source/sink setup."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import usgs_series


class UsgsSeriesTests(unittest.TestCase):
    def test_station_id_normalization_preserves_leading_zeros(self):
        self.assertEqual(
            usgs_series._normalize_usgs_station_id("1358000"),
            "01358000",
        )
        self.assertEqual(
            usgs_series._normalize_usgs_station_id("01358000.0"),
            "01358000",
        )

    def test_hudson_chunks_are_nonoverlapping_and_cover_period(self):
        chunks = list(
            usgs_series._iter_hudson_chunks(
                "2020-01-01",
                "2020-01-07",
                chunk_days=3,
            )
        )

        self.assertEqual(
            [(start.day, end.day) for _, start, end in chunks],
            [(1, 3), (4, 6), (7, 7)],
        )

    def test_merge_usgs_series_sorts_and_deduplicates(self):
        first = pd.Series(
            [1.0, 2.0],
            index=pd.to_datetime(
                ["2020-01-02T00:00Z", "2020-01-03T00:00Z"]
            ),
        )
        second = pd.Series(
            [4.0, 3.0],
            index=pd.to_datetime(
                ["2020-01-03T00:00Z", "2020-01-01T00:00Z"]
            ),
        )

        merged = usgs_series._merge_usgs_series([first, second])

        np.testing.assert_array_equal(merged.to_numpy(), [3.0, 1.0, 2.0])
        self.assertTrue(merged.index.is_monotonic_increasing)

    def test_interpolation_uses_short_gaps_and_preserves_other_values(self):
        observations = pd.Series(
            [10.0, 14.0, 30.0],
            index=pd.to_datetime(
                [
                    "2020-01-01T01:00Z",
                    "2020-01-01T05:00Z",
                    "2020-01-01T20:00Z",
                ]
            ),
        )
        target_time = pd.date_range(
            "2020-01-01T00:00Z",
            periods=22,
            freq="h",
        )
        original = np.full(len(target_time), -1.0)

        blended, use_usgs = usgs_series._interpolate_usgs_with_original_fallback(
            series=observations,
            target_time=target_time,
            original_values=original,
        )

        self.assertEqual(blended[0], -1.0)
        self.assertEqual(blended[1], 10.0)
        self.assertEqual(blended[3], 12.0)
        self.assertEqual(blended[10], -1.0)
        self.assertEqual(blended[20], 30.0)
        self.assertEqual(blended[21], -1.0)
        self.assertTrue(np.all(use_usgs[1:6]))
        self.assertFalse(use_usgs[10])

    def test_hudson_adaptive_download_splits_failed_window(self):
        successful_piece = pd.Series(
            [1.0],
            index=pd.to_datetime(["2020-01-01T00:00Z"]),
        )

        def fake_download(
            station_id,
            parameter_id,
            window_start,
            window_end,
            cache_dir,
        ):
            del station_id, parameter_id, cache_dir
            span_days = (window_end - window_start).days + 1
            return successful_piece if span_days <= 5 else None

        with patch.object(
            usgs_series,
            "_download_hudson_window",
            side_effect=fake_download,
        ) as downloader:
            pieces = usgs_series._download_hudson_window_adaptive(
                station_id="01358000",
                parameter_id="00060",
                window_start=pd.Timestamp("2020-01-01", tz="UTC"),
                window_end=pd.Timestamp("2020-01-10", tz="UTC"),
                cache_dir=None,
                retry_chunk_days=20,
                min_chunk_days=5,
            )

        self.assertEqual(len(pieces), 2)
        self.assertEqual(downloader.call_count, 3)


if __name__ == "__main__":
    unittest.main()
