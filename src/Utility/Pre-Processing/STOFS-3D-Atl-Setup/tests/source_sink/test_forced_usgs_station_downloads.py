"""Tests for strict retries of manually linked USGS stations."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import download_usgs
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    FORCED_STATION_MAX_MISSING_GAP_DAYS,
    FORCED_STATION_RETRY_POLICY,
    OPTIONAL_FLOW_STATIONS,
    OPTIONAL_TEMPERATURE_STATIONS,
    REQUIRED_FLOW_STATIONS,
    REQUIRED_TEMPERATURE_STATIONS,
)


class ForcedUsgsStationDownloadTests(unittest.TestCase):
    def test_forced_station_retry_policy_has_adaptive_windows(self):
        self.assertEqual(
            FORCED_STATION_RETRY_POLICY["primary_chunk_days"], 100
        )
        self.assertEqual(
            FORCED_STATION_RETRY_POLICY["retry_chunk_days"],
            [20, 10, 5],
        )
        self.assertEqual(FORCED_STATION_MAX_MISSING_GAP_DAYS, 7)

    def _download_with_missing_primary_station(
        self,
        retry_result,
        period_days=2,
        fetch_results=None,
        observations=None,
        requirement="required",
    ):
        station_id = "01358000"
        if observations is None:
            observations = pd.Series(
                [4.0, 5.0],
                index=pd.to_datetime(
                    ["2020-01-01T00:00Z", "2020-01-02T00:00Z"]
                ),
            )
        metadata = pd.DataFrame(
            {
                "site_no": [station_id],
                "station_nm": ["Hudson River at Albany"],
                "dec_long_va": [-73.75],
                "dec_lat_va": [42.65],
            }
        )
        if fetch_results is None:
            fetch_results = [pd.DataFrame(), pd.DataFrame()]

        with tempfile.TemporaryDirectory() as directory:
            with (
                patch.object(download_usgs, "make_session"),
                patch.object(download_usgs, "usgs_site_info", return_value=metadata),
                patch.object(
                    download_usgs,
                    "fetch_iv_rdb",
                    side_effect=fetch_results,
                ),
                patch(
                    "stofs3d_setup.ops.Source_sink.Replace_with_USGS."
                    "usgs_series.download_usgs_series",
                    return_value=retry_result(observations),
                ) as adaptive_retry,
                patch.object(download_usgs, "write_manifest"),
                patch.object(download_usgs, "save_cache_parquet"),
            ):
                records = download_usgs.download_stations(
                    param_id="00060",
                    station_ids=[station_id],
                    cache_fname=str(Path(directory) / "flow.pq"),
                    datelist=pd.date_range(
                        "2020-01-01", periods=period_days
                    ),
                    required_station_ids=(
                        [station_id] if requirement == "required" else []
                    ),
                    optional_station_ids=(
                        [station_id] if requirement == "optional" else []
                    ),
                    required_station_retry_policy=FORCED_STATION_RETRY_POLICY,
                    required_station_max_missing_gap_days=7,
                    required_station_start="2020-01-01",
                    required_station_end=(
                        pd.Timestamp("2020-01-01")
                        + pd.Timedelta(days=period_days - 1)
                    ),
                    warning_log_path=Path(directory) / "warnings.log",
                    warning_input_files=(
                        "Source_sink/vsource.th",
                        "Source_sink/source.nc",
                    ),
                    warning_forcing_name="streamflow",
                )
                warning_log = (Path(directory) / "warnings.log").read_text(
                    encoding="utf-8"
                ) if (Path(directory) / "warnings.log").exists() else ""

        return records, adaptive_retry, observations, warning_log

    def test_required_station_uses_adaptive_second_fallback(self):
        records, adaptive_retry, observations, _ = (
            self._download_with_missing_primary_station(
                lambda series: series
            )
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].station_info["id"], "01358000")
        self.assertEqual(adaptive_retry.call_count, 1)
        self.assertEqual(
            adaptive_retry.call_args.kwargs["policy"],
            FORCED_STATION_RETRY_POLICY,
        )
        pd.testing.assert_series_equal(
            pd.Series(
                records[0].df["value"].to_numpy(),
                index=pd.DatetimeIndex(records[0].df["date"]),
            ),
            observations,
            check_freq=False,
            check_names=False,
        )

    def test_required_station_raises_after_both_fallbacks_fail(self):
        with self.assertRaisesRegex(
            download_usgs.RequiredStationDownloadError,
            "after both fallback methods.*01358000",
        ):
            self._download_with_missing_primary_station(lambda series: None)

    def test_gappy_required_station_gets_both_fallbacks(self):
        station_id = "01358000"
        gappy_data = pd.DataFrame(
            {
                "site_no": [station_id, station_id],
                "date": pd.to_datetime(["2020-01-01", "2020-01-10"]),
                "value": [4.0, 5.0],
                "tz_cd": ["UTC", "UTC"],
            }
        )
        complete_observations = pd.Series(
            [4.0, 4.5, 5.0],
            index=pd.to_datetime(
                [
                    "2020-01-01T00:00Z",
                    "2020-01-05T00:00Z",
                    "2020-01-10T00:00Z",
                ]
            ),
        )

        records, adaptive_retry, _, _ = self._download_with_missing_primary_station(
            lambda series: complete_observations,
            period_days=10,
            fetch_results=[gappy_data, gappy_data],
            observations=complete_observations,
        )

        self.assertEqual(adaptive_retry.call_count, 1)
        self.assertEqual(records[0].station_info["id"], station_id)
        self.assertEqual(len(records[0].df), 3)

    def test_station_still_over_gap_threshold_raises_after_retry(self):
        station_id = "01358000"
        gappy_data = pd.DataFrame(
            {
                "site_no": [station_id, station_id],
                "date": pd.to_datetime(["2020-01-01", "2020-01-10"]),
                "value": [4.0, 5.0],
                "tz_cd": ["UTC", "UTC"],
            }
        )
        gappy_observations = pd.Series(
            [4.0, 5.0],
            index=pd.to_datetime(
                ["2020-01-01T00:00Z", "2020-01-10T00:00Z"]
            ),
        )

        with self.assertRaisesRegex(
            download_usgs.RequiredStationDownloadError,
            "largest missing duration.*exceeds",
        ):
            self._download_with_missing_primary_station(
                lambda series: gappy_observations,
                period_days=10,
                fetch_results=[gappy_data, gappy_data],
                observations=gappy_observations,
            )

    def test_optional_station_logs_warning_without_raising(self):
        records, adaptive_retry, _, warning_log = (
            self._download_with_missing_primary_station(
                lambda series: None,
                requirement="optional",
            )
        )

        self.assertEqual(records, [])
        self.assertEqual(adaptive_retry.call_count, 1)
        self.assertIn("Forced USGS streamflow replacement", warning_log)
        self.assertIn("01358000", warning_log)
        self.assertIn("Source_sink/vsource.th", warning_log)
        self.assertIn("Source_sink/source.nc", warning_log)

    def test_manual_soft_link_station_requirements(self):
        self.assertIn("07375175", OPTIONAL_FLOW_STATIONS)
        self.assertNotIn("07375175", REQUIRED_FLOW_STATIONS)
        self.assertIn("02244440", OPTIONAL_TEMPERATURE_STATIONS)
        self.assertNotIn("02244440", REQUIRED_TEMPERATURE_STATIONS)
        self.assertIn("07381490", REQUIRED_FLOW_STATIONS)


if __name__ == "__main__":
    unittest.main()
