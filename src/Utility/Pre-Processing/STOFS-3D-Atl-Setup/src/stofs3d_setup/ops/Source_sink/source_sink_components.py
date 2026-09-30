"""Aligned array operations for SCHISM source/sink forcing."""

import numpy as np

from pylib_experimental.schism_file import source_sink, TimeHistory


def _time_array(th: TimeHistory) -> np.ndarray:
    """Return TimeHistory time values as a one-dimensional array."""
    return np.asarray(th.time).reshape(-1)


def _data_array(th: TimeHistory) -> np.ndarray:
    """Return TimeHistory values as a two-dimensional float array."""
    data = np.asarray(th.data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    return data


def _build_timehistory(
    time: np.ndarray,
    data: np.ndarray,
    element_ids: list[int],
) -> TimeHistory | None:
    """Create a TimeHistory object, or return None when no elements remain."""
    element_ids = [int(ele) for ele in element_ids]

    if len(element_ids) == 0:
        return None

    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if data.shape[1] != len(element_ids):
        raise ValueError(
            "TimeHistory data/element mismatch: "
            f"{data.shape[1]} columns for {len(element_ids)} elements"
        )

    return TimeHistory(
        data_array=np.c_[np.asarray(time).reshape(-1), data],
        columns=[str(ele) for ele in element_ids],
    )


def _copy_source_components(base_ss):
    """Extract source element IDs, vsource data, and msource data."""
    source_eles = [int(ele) for ele in np.asarray(base_ss.source_eles).reshape(-1)]

    if base_ss.vsource is None:
        if source_eles:
            raise ValueError("base_ss has source elements but no vsource")
        return [], None, []

    vsource_time = _time_array(base_ss.vsource)
    vsource_data = _data_array(base_ss.vsource).copy()

    if vsource_data.shape[1] != len(source_eles):
        raise ValueError(
            f"base_ss source count is {len(source_eles)}, "
            f"but vsource has {vsource_data.shape[1]} columns"
        )

    msource_data_list = []
    for msource in base_ss.msource or []:
        msource_data = _data_array(msource).copy()
        if msource_data.shape[1] != len(source_eles):
            raise ValueError(
                f"base_ss source count is {len(source_eles)}, "
                f"but an msource tracer has {msource_data.shape[1]} columns"
            )
        msource_data_list.append(
            (_time_array(msource), msource_data)
        )

    return source_eles, (vsource_time, vsource_data), msource_data_list


def _copy_sink_components(base_ss):
    """Extract sink element IDs and vsink data."""
    sink_eles_raw = getattr(base_ss, "sink_eles", None)

    if sink_eles_raw is None:
        sink_eles = []
    else:
        sink_eles = [int(ele) for ele in np.asarray(sink_eles_raw).reshape(-1)]

    if base_ss.vsink is None:
        if sink_eles:
            raise ValueError("base_ss has sink elements but no vsink")
        return [], None

    vsink_time = _time_array(base_ss.vsink)
    vsink_data = _data_array(base_ss.vsink).copy()

    if vsink_data.shape[1] != len(sink_eles):
        raise ValueError(
            f"base_ss sink count is {len(sink_eles)}, "
            f"but vsink has {vsink_data.shape[1]} columns"
        )

    return sink_eles, (vsink_time, vsink_data)


def _move_source_element(
    source_eles: list[int],
    old_ele: int,
    target_ele: int,
) -> None:
    """Relabel one source element without changing its time series."""
    if old_ele not in source_eles:
        raise ValueError(f"Source element {old_ele} is not in base_ss")

    if target_ele != old_ele and target_ele in source_eles:
        raise ValueError(
            f"Cannot move source element {old_ele} to {target_ele}: "
            "the target element already contains a source"
        )

    source_eles[source_eles.index(old_ele)] = target_ele


def _move_sink_element(
    sink_eles: list[int],
    old_ele: int,
    target_ele: int,
) -> None:
    """Relabel one sink element without changing its time series."""
    if old_ele not in sink_eles:
        raise ValueError(f"Sink element {old_ele} is not in base_ss")

    if target_ele != old_ele and target_ele in sink_eles:
        raise ValueError(
            f"Cannot move sink element {old_ele} to {target_ele}: "
            "the target element already contains a sink"
        )

    sink_eles[sink_eles.index(old_ele)] = target_ele


def _add_negative_source_part_to_sink(
    target_ele: int,
    source_column_idx: int,
    source_data: np.ndarray,
    source_time: np.ndarray,
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
) -> tuple[list[int], tuple[np.ndarray, np.ndarray] | None, int, float]:
    """
    Move negative source values to a sink at the same element.

    Returns updated sink information, number of negative records, and minimum
    negative flow.
    """
    raw_source = source_data[:, source_column_idx].copy()
    negative_part = np.minimum(raw_source, 0.0)
    n_negative = int(np.sum(negative_part < 0.0))
    min_negative = float(np.min(negative_part)) if len(negative_part) else 0.0

    source_data[:, source_column_idx] = np.maximum(raw_source, 0.0)

    if n_negative == 0:
        return sink_eles, sink_time_and_data, 0, min_negative

    if sink_time_and_data is None:
        sink_time = source_time.copy()
        sink_data = negative_part.reshape(-1, 1)
        sink_eles = [int(target_ele)]
        return sink_eles, (sink_time, sink_data), n_negative, min_negative

    sink_time, sink_data = sink_time_and_data

    if len(sink_time) != len(source_time) or not np.allclose(
        np.asarray(sink_time, dtype=float),
        np.asarray(source_time, dtype=float),
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            "Cannot add artificial-island negative sink because vsource and "
            "vsink time arrays differ"
        )

    if target_ele in sink_eles:
        sink_idx = sink_eles.index(target_ele)
        sink_data[:, sink_idx] += negative_part
    else:
        sink_eles.append(int(target_ele))
        sink_data = np.column_stack((sink_data, negative_part))

    return sink_eles, (sink_time, sink_data), n_negative, min_negative


def _remove_source_columns(
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    msource_data_list: list[tuple[np.ndarray, np.ndarray]],
    remove_eles: set[int],
) -> tuple[
    list[int],
    tuple[np.ndarray, np.ndarray] | None,
    list[tuple[np.ndarray, np.ndarray]],
]:
    """Remove selected source elements from vsource and all msource tracers."""
    if not remove_eles or source_time_and_data is None:
        return source_eles, source_time_and_data, msource_data_list

    keep = np.array(
        [ele not in remove_eles for ele in source_eles],
        dtype=bool,
    )

    source_time, source_data = source_time_and_data
    source_eles = [ele for ele, keep_it in zip(source_eles, keep) if keep_it]
    source_data = source_data[:, keep]

    updated_msource = []
    for tracer_time, tracer_data in msource_data_list:
        updated_msource.append((tracer_time, tracer_data[:, keep]))

    if len(source_eles) == 0:
        return [], None, []

    return (
        source_eles,
        (source_time, source_data),
        updated_msource,
    )


def _remove_sink_columns(
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    remove_eles: set[int],
) -> tuple[list[int], tuple[np.ndarray, np.ndarray] | None]:
    """Remove selected sink elements from vsink."""
    if not remove_eles or sink_time_and_data is None:
        return sink_eles, sink_time_and_data

    keep = np.array(
        [ele not in remove_eles for ele in sink_eles],
        dtype=bool,
    )

    sink_time, sink_data = sink_time_and_data
    sink_eles = [ele for ele, keep_it in zip(sink_eles, keep) if keep_it]
    sink_data = sink_data[:, keep]

    if len(sink_eles) == 0:
        return [], None

    return sink_eles, (sink_time, sink_data)


def _build_source_sink(
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    msource_data_list: list[tuple[np.ndarray, np.ndarray]],
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
) -> source_sink:
    """Build a new source_sink object from patched arrays."""
    if source_time_and_data is None:
        vsource = None
        msource_list = None
    else:
        source_time, source_data = source_time_and_data
        vsource = _build_timehistory(source_time, source_data, source_eles)

        msource_list = []
        for tracer_time, tracer_data in msource_data_list:
            msource_list.append(
                _build_timehistory(tracer_time, tracer_data, source_eles)
            )

    if sink_time_and_data is None:
        vsink = None
    else:
        sink_time, sink_data = sink_time_and_data
        vsink = _build_timehistory(sink_time, sink_data, sink_eles)

    return source_sink(
        vsource=vsource,
        vsink=vsink,
        msource=msource_list,
    )

def _check_matching_time(
    reference_time: np.ndarray,
    candidate_time: np.ndarray,
    label: str,
) -> None:
    """Require two SCHISM TimeHistory time arrays to match."""
    reference_time = np.asarray(reference_time, dtype=float).reshape(-1)
    candidate_time = np.asarray(candidate_time, dtype=float).reshape(-1)

    if len(reference_time) != len(candidate_time) or not np.allclose(
        reference_time,
        candidate_time,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            f"Time array mismatch while restoring {label}: "
            f"base={len(reference_time)} records, "
            f"original={len(candidate_time)} records"
        )

def _append_source_column(
    source_eles: list[int],
    source_time_and_data,
    msource_data_list,
    target_ele: int,
    source_time: np.ndarray,
    flow: np.ndarray,
    tracer_columns: list[np.ndarray],
):
    """Append one new source and corresponding msource tracer columns."""
    if target_ele in source_eles:
        raise ValueError(
            f"Target element {target_ele} already contains a source"
        )

    flow = np.asarray(flow, dtype=float).reshape(-1)
    source_time = np.asarray(source_time, dtype=float).reshape(-1)

    if len(flow) != len(source_time):
        raise ValueError("New source flow length does not match model time")

    if source_time_and_data is None:
        combined_source_time = source_time.copy()
        combined_source_data = flow.reshape(-1, 1)
    else:
        combined_source_time, combined_source_data = source_time_and_data
        _check_matching_time(
            combined_source_time,
            source_time,
            f"new source element {target_ele}",
        )
        combined_source_data = np.column_stack(
            (combined_source_data, flow)
        )

    if len(msource_data_list) != len(tracer_columns):
        raise ValueError(
            "New source tracer count does not match base msource tracer count"
        )

    updated_msource = []
    for tracer_idx, (
        (tracer_time, tracer_data),
        tracer_column,
    ) in enumerate(zip(msource_data_list, tracer_columns)):
        _check_matching_time(
            tracer_time,
            source_time,
            f"new source msource tracer {tracer_idx + 1}",
        )
        tracer_column = np.asarray(tracer_column, dtype=float).reshape(-1)
        if len(tracer_column) != len(source_time):
            raise ValueError(
                f"Tracer {tracer_idx + 1} length does not match model time"
            )
        updated_msource.append(
            (
                tracer_time,
                np.column_stack((tracer_data, tracer_column)),
            )
        )

    source_eles.append(int(target_ele))
    return (
        source_eles,
        (combined_source_time, combined_source_data),
        updated_msource,
    )


def _append_sink_column(
    sink_eles: list[int],
    sink_time_and_data,
    target_ele: int,
    sink_time: np.ndarray,
    sink_values: np.ndarray,
):
    """Append one new sink column."""
    if target_ele in sink_eles:
        raise ValueError(
            f"Target element {target_ele} already contains a sink"
        )

    sink_time = np.asarray(sink_time, dtype=float).reshape(-1)
    sink_values = np.asarray(sink_values, dtype=float).reshape(-1)

    if len(sink_time) != len(sink_values):
        raise ValueError("New sink length does not match model time")

    if sink_time_and_data is None:
        combined_time = sink_time.copy()
        combined_data = sink_values.reshape(-1, 1)
    else:
        combined_time, combined_data = sink_time_and_data
        _check_matching_time(
            combined_time,
            sink_time,
            f"new sink element {target_ele}",
        )
        combined_data = np.column_stack(
            (combined_data, sink_values)
        )

    sink_eles.append(int(target_ele))
    return sink_eles, (combined_time, combined_data)
