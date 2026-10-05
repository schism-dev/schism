"""Characterization tests for automatic source-temperature support."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS import source_temperature
from stofs3d_setup.ops.Source_sink.Replace_with_USGS import hard_reset_config
from stofs3d_setup.ops.Source_sink.Replace_with_USGS import (
    nwm_usgs_station_search,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    NWM_TO_USGS_FLOW_ADJUSTMENT,
    NWM_TO_USGS_TEMPERATURE_SEARCH,
    TEMPERATURE_SEARCH_DIRECTIONS,
)
from stofs3d_setup.ops.Source_sink.correction_config import (
    load_source_sink_corrections,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
)


class _Grid:
    xctr = np.array([0.0])
    yctr = np.array([0.0])

    def compute_ctr(self):
        pass


class SourceTemperatureTests(unittest.TestCase):
    def test_shared_station_search_applies_workflow_specific_links(self):
        network = pd.DataFrame(
            {
                "featureID": [101, 202],
                "gages": [None, None],
            }
        )
        with TemporaryDirectory() as cache_dir:
            with (
                patch.object(
                    nwm_usgs_station_search,
                    "prepare_usgs_stations",
                    return_value=(
                        np.array(["01234567"]),
                        np.array([[-75.0, 39.0]]),
                    ),
                ),
                patch.object(
                    nwm_usgs_station_search,
                    "preprocess_nwm_shp",
                    return_value=network.copy(),
                ),
                patch.object(
                    nwm_usgs_station_search,
                    "associate_poi_with_nwm",
                    side_effect=lambda frame, **kwargs: frame,
                ),
            ):
                result = nwm_usgs_station_search.prepare_nwm_usgs_station_search(
                    states=["DE"],
                    nwm_shapefile="unused.shp",
                    usgs_cache_folder=cache_dir,
                    manual_feature_station_links={202: "01359139"},
                    invalid_station_diagnostic_path=None,
                )

        self.assertEqual(result.loc[1, "gages"], "01359139")
        self.assertIsNone(result.loc[0, "gages"])

    def test_source_mapping_loader_normalizes_elements_and_feature_ids(self):
        with TemporaryDirectory() as directory:
            mapping_file = Path(directory) / "sources.json"
            mapping_file.write_text(
                json.dumps({"10": [1, "2"], "20": 3}),
                encoding="utf-8",
            )

            mapping = source_temperature._load_source_feature_mapping(
                directory
            )

        self.assertEqual(mapping, {10: [1, 2], 20: [3]})

    def test_source_mapping_loader_requires_sources_json(self):
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "sources.json"):
                source_temperature._load_source_feature_mapping(directory)

    def test_explicit_hudson_station_mappings_are_parameter_specific(self):
        hard_reset_file = Path(hard_reset_config.__file__).with_name(
            "hard_reset_sources_with_usgs.yml"
        )
        corrections = load_source_sink_corrections(hard_reset_file)
        hard_reset_points = hard_reset_config.source_override_points(
            corrections
        )
        hudson = next(
            point
            for point in hard_reset_points
            if point["name"] == "Hudson River"
        )

        self.assertEqual(hudson["flow_station_id"], "01358000")
        self.assertEqual(hudson["temperature_station_id"], "01359139")
        self.assertEqual(NWM_TO_USGS_TEMPERATURE_SEARCH[6186156], "01359139")
        self.assertEqual(NWM_TO_USGS_FLOW_ADJUSTMENT[6186156], "01358000")
        self.assertNotIn(9643251, NWM_TO_USGS_TEMPERATURE_SEARCH)

    def test_manual_feature_station_links_preserve_both_scopes(self):
        flow_links = {
            19406836: "07381490",
            15708755: "02489500",
            18928090: "07375175",
            19269176: "07374000",
            9643251: "02172035",
            16665157: "02244040",
            16665419: "02244440",
            2590217: "01463500",
            6186156: "01358000",
        }
        temperature_additions = {
            16665157: "02244040",
            16665419: "02244440",
            2590217: "01463500",
            6186156: "01359139",
        }

        self.assertEqual(NWM_TO_USGS_FLOW_ADJUSTMENT, flow_links)
        self.assertEqual(
            NWM_TO_USGS_TEMPERATURE_SEARCH,
            temperature_additions,
        )

    def test_station_search_preserves_station_feature_pair(self):
        def add_candidate(**kwargs):
            kwargs["vsource"].usgs_st.append(
                SimpleNamespace(st_id="01234567", nearby_nwm_fid=202)
            )

        with patch.object(
            source_temperature,
            "find_usgs_along_nwm",
            side_effect=add_candidate,
        ):
            find_candidates = (
                source_temperature._find_temperature_station_candidates
            )
            candidates = find_candidates(
                source_eles=[1],
                source_element_to_fids={1: [101]},
                source_time_and_data=(
                    np.array([0.0]),
                    np.array([[1.0]]),
                ),
                xctr=np.array([0.0]),
                yctr=np.array([0.0]),
                nwm_shp=None,
            )

        self.assertEqual(
            candidates[1],
            [
                source_temperature.TemperatureStationCandidate(
                    station_id="01234567",
                    nwm_feature_id=202,
                )
            ],
        )
        self.assertEqual(
            TEMPERATURE_SEARCH_DIRECTIONS,
            ("upstream", "downstream"),
        )

    def test_temperature_search_can_include_downstream_without_flow_change(self):
        searched_directions = []

        def add_candidate(**kwargs):
            searched_directions.append(kwargs["iupstream"])
            kwargs["vsource"].usgs_st.append(
                SimpleNamespace(
                    st_id=("upstream" if kwargs["iupstream"] else "downstream"),
                    nearby_nwm_fid=202 if kwargs["iupstream"] else 303,
                )
            )

        with patch.object(
            source_temperature,
            "find_usgs_along_nwm",
            side_effect=add_candidate,
        ):
            candidates = source_temperature._find_temperature_station_candidates(
                source_eles=[1],
                source_element_to_fids={1: [101]},
                source_time_and_data=(
                    np.array([0.0]),
                    np.array([[1.0]]),
                ),
                xctr=np.array([0.0]),
                yctr=np.array([0.0]),
                nwm_shp=None,
            )

        self.assertEqual(searched_directions, [True, False])
        self.assertEqual(
            [candidate.station_id for candidate in candidates[1]],
            ["upstream", "downstream"],
        )

    def test_manual_temperature_station_is_first_usable_candidate(self):
        def add_candidates(**kwargs):
            kwargs["vsource"].usgs_st.extend(
                [
                    SimpleNamespace(st_id="automatic", nearby_nwm_fid=202),
                    SimpleNamespace(st_id="manual", nearby_nwm_fid=101),
                ]
            )

        target_time = pd.date_range(
            "2020-01-01", periods=2, freq="1h", tz="UTC"
        )
        with (
            patch.dict(
                source_temperature.NWM_TO_USGS_TEMPERATURE_SEARCH,
                {101: "manual"},
            ),
            patch.object(
                source_temperature,
                "find_usgs_along_nwm",
                side_effect=add_candidates,
            ),
        ):
            candidates = (
                source_temperature._find_temperature_station_candidates(
                    source_eles=[1],
                    source_element_to_fids={1: [101]},
                    source_time_and_data=(
                        np.array([0.0]),
                        np.array([[1.0]]),
                    ),
                    xctr=np.array([0.0]),
                    yctr=np.array([0.0]),
                    nwm_shp=None,
                )[1]
            )

        replacement, station_id, _ = (
            source_temperature._first_usable_temperature(
                candidates=candidates,
                temperature_by_station={
                    "automatic": pd.Series([5.0, 5.0], index=target_time),
                    "manual": pd.Series([15.0, 16.0], index=target_time),
                },
                target_time=target_time,
                original_values=np.array([-9999.0, -9999.0]),
            )
        )

        self.assertEqual(station_id, "manual")
        np.testing.assert_array_equal(replacement, [15.0, 16.0])

    def test_public_stage_changes_temperature_without_changing_flow(self):
        time = np.array([0.0, 3600.0])
        original = _build_source_sink(
            source_eles=[1],
            source_time_and_data=(time, np.array([[2.0], [3.0]])),
            msource_data_list=[
                (time, np.array([[10.0], [11.0]])),
                (time, np.array([[0.0], [0.0]])),
            ],
            sink_eles=[],
            sink_time_and_data=None,
        )

        def replace_temperature(**kwargs):
            tracers = kwargs["msource_data_list"]
            tracers[0][1][:, 0] = [20.0, 21.0]
            return tracers, 1

        with patch.object(
            source_temperature,
            "_apply_automatic_temperature_replacements",
            side_effect=replace_temperature,
        ):
            corrected, count = (
                source_temperature.replace_source_temperatures_with_usgs(
                    base_ss=original,
                    hgrid=_Grid(),
                    source_mapping_dir="unused",
                    start_time="2020-01-01",
                    usgs_cache_folder="unused",
                    nwm_shapefile="unused",
                )
            )

        self.assertEqual(count, 1)
        self.assertEqual(np.asarray(corrected.source_eles).tolist(), [1])
        np.testing.assert_array_equal(
            corrected.vsource.df.values,
            [[2.0], [3.0]],
        )
        np.testing.assert_array_equal(
            corrected.msource[0].df.values,
            [[20.0], [21.0]],
        )
        np.testing.assert_array_equal(
            original.msource[0].df.values,
            [[10.0], [11.0]],
        )

    def test_discharge_weighted_temperature_uses_associated_feature_flow(self):
        target_time = pd.date_range(
            "2020-01-01", periods=3, freq="1h", tz="UTC"
        )
        candidates = [
            source_temperature.TemperatureStationCandidate("A", 101),
            source_temperature.TemperatureStationCandidate("B", 202),
        ]
        temperatures = {
            "A": pd.Series([10.0, 10.0, 10.0], index=target_time),
            "B": pd.Series([20.0, 20.0, 20.0], index=target_time),
        }
        discharges = {
            101: np.array([1.0, 3.0, 1.0]),
            202: np.array([3.0, 1.0, 1.0]),
        }

        result, diagnostics = (
            source_temperature._discharge_weighted_temperature(
                candidates=candidates,
                temperature_by_station=temperatures,
                discharge_by_feature=discharges,
                target_time=target_time,
            )
        )

        np.testing.assert_allclose(result, [17.5, 12.5, 15.0])
        self.assertEqual(diagnostics["complete"], 1)

    def test_discharge_weighted_temperature_has_no_spatial_threshold(self):
        target_time = pd.date_range(
            "2020-01-01", periods=2, freq="1h", tz="UTC"
        )
        candidates = [
            source_temperature.TemperatureStationCandidate("creek", 1),
            source_temperature.TemperatureStationCandidate("main", 2),
        ]
        temperatures = {
            "creek": pd.Series([8.0, 9.0], index=target_time),
        }
        discharges = {
            1: np.array([1.0, 1.0]),
            2: np.array([99.0, 99.0]),
        }

        result, diagnostics = (
            source_temperature._discharge_weighted_temperature(
                candidates=candidates,
                temperature_by_station=temperatures,
                discharge_by_feature=discharges,
                target_time=target_time,
            )
        )

        np.testing.assert_array_equal(result, [8.0, 9.0])
        self.assertAlmostEqual(
            diagnostics["minimum_discharge_coverage"], 0.01
        )

    def test_incomplete_weighted_temperature_preserves_ambient_column(self):
        target_time = pd.date_range(
            "2020-01-01", periods=3, freq="1h", tz="UTC"
        )
        candidates = [
            source_temperature.TemperatureStationCandidate("A", 101),
        ]
        temperatures = {
            "A": pd.Series([10.0], index=target_time[:1]),
        }

        result, diagnostics = (
            source_temperature._discharge_weighted_temperature(
                candidates=candidates,
                temperature_by_station=temperatures,
                discharge_by_feature={101: np.ones(3)},
                target_time=target_time,
            )
        )

        self.assertIsNone(result)
        self.assertEqual(diagnostics["complete"], 0)

    def test_atomic_temperature_blend_rejects_partial_station_series(self):
        target_time = pd.date_range(
            "2020-01-01", periods=3, freq="1h", tz="UTC"
        )
        original = np.full(3, -9999.0)
        partial = pd.Series([10.0], index=target_time[:1])

        blended, use_usgs, applied = (
            source_temperature._blend_complete_temperature_or_preserve(
                series=partial,
                target_time=target_time,
                original_values=original,
            )
        )

        np.testing.assert_array_equal(blended, original)
        np.testing.assert_array_equal(use_usgs, [False, False, False])
        self.assertFalse(applied)


if __name__ == "__main__":
    unittest.main()
