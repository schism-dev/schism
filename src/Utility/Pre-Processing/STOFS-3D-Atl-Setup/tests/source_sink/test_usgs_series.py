"""Characterization tests for USGS series handling used by source/sink setup."""

import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from stofs3d_setup.ops.USGS.ogc import series as usgs_series


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

    def test_date_windows_are_nonoverlapping_and_cover_period(self):
        chunks = list(
            usgs_series._iter_date_windows(
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

    def test_adaptive_download_uses_configured_smaller_windows(self):
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
            cache_prefix,
        ):
            del station_id, parameter_id, cache_dir, cache_prefix
            span_days = (window_end - window_start).days + 1
            return usgs_series.UsgsWindowResult(
                series=successful_piece if span_days <= 5 else None,
                retryable=span_days > 5,
            )

        with patch.object(
            usgs_series,
            "_download_usgs_window",
            side_effect=fake_download,
        ) as downloader:
            pieces = usgs_series._download_usgs_window_adaptive(
                station_id="01358000",
                parameter_id="00060",
                window_start=pd.Timestamp("2020-01-01", tz="UTC"),
                window_end=pd.Timestamp("2020-01-10", tz="UTC"),
                cache_dir=None,
                cache_prefix="test_usgs",
                retry_chunk_days=(20, 10, 5),
            )

        self.assertEqual(len(pieces), 2)
        self.assertEqual(downloader.call_count, 3)

    def test_successful_empty_window_is_not_retried(self):
        with patch.object(
            usgs_series,
            "_download_usgs_window",
            return_value=usgs_series.UsgsWindowResult(
                series=None,
                retryable=False,
            ),
        ) as downloader:
            pieces = usgs_series._download_usgs_window_adaptive(
                station_id="02244440",
                parameter_id="00010",
                window_start=pd.Timestamp("2020-01-01", tz="UTC"),
                window_end=pd.Timestamp("2020-04-09", tz="UTC"),
                cache_dir=None,
                cache_prefix="test_usgs",
                retry_chunk_days=(20, 10, 5),
            )

        self.assertEqual(pieces, [])
        downloader.assert_called_once()

    def test_observation_free_period_only_uses_primary_chunks(self):
        empty_window = usgs_series.UsgsWindowResult(
            series=None,
            retryable=False,
        )
        policy = usgs_series.UsgsDownloadPolicy(
            primary_chunk_days=3,
            retry_chunk_days=(2, 1),
            cache_prefix="test_usgs",
        )

        with TemporaryDirectory() as cache_dir:
            with patch.object(
                usgs_series,
                "_download_usgs_window",
                return_value=empty_window,
            ) as downloader:
                result = usgs_series.download_usgs_series(
                    station_id="02244440",
                    parameter_id="00010",
                    start_time="2020-01-01",
                    end_time="2020-01-10",
                    usgs_cache_folder=cache_dir,
                    policy=policy,
                )

        self.assertIsNone(result)
        # The downloader pads the requested period by one day at both ends:
        # 12 days / 3-day primary chunks, with no adaptive subdivisions.
        self.assertEqual(downloader.call_count, 4)

    def test_download_policy_requires_decreasing_positive_windows(self):
        with self.assertRaisesRegex(ValueError, "strictly decreasing"):
            usgs_series.UsgsDownloadPolicy(
                primary_chunk_days=20,
                retry_chunk_days=(20, 5),
            )
        with self.assertRaisesRegex(ValueError, "positive"):
            usgs_series.UsgsDownloadPolicy(retry_chunk_days=(0,))

    def test_default_download_policy_is_loaded_from_shared_yaml(self):
        self.assertEqual(
            usgs_series.DEFAULT_USGS_DOWNLOAD_POLICY.primary_chunk_days,
            100,
        )
        self.assertEqual(
            usgs_series.DEFAULT_USGS_DOWNLOAD_POLICY.retry_chunk_days,
            (20, 10, 5),
        )

    def test_named_series_uses_configured_policy_without_name_branch(self):
        policy = {
            "primary_chunk_days": 100,
            "retry_chunk_days": [20, 10, 5],
            "cache_prefix": "custom_usgs",
        }
        downloaded = pd.Series(
            [1.0],
            index=pd.to_datetime(["2020-01-01T00:00Z"]),
        )

        with patch.object(
            usgs_series,
            "download_usgs_series",
            return_value=downloaded,
        ) as downloader:
            result, station_id = usgs_series._download_configured_usgs_series(
                station_id="1358000",
                parameter_id="00060",
                start_time="2020-01-01",
                model_time=np.array([0.0, 3600.0]),
                usgs_cache_folder="unused",
                download_policy=policy,
            )

        self.assertIs(result, downloaded)
        self.assertEqual(station_id, "01358000")
        self.assertEqual(downloader.call_args.kwargs["policy"], policy)

    def test_download_policy_rejects_unknown_configuration(self):
        with self.assertRaisesRegex(ValueError, "unsupported_setting"):
            usgs_series._as_usgs_download_policy(
                {"unsupported_setting": 10}
            )


if __name__ == "__main__":
    unittest.main()
