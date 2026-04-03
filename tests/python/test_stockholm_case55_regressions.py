from __future__ import annotations

import re
from pathlib import Path

import pytest

import dtcc_mesher as dm


_CASE55_BAD_DOMAINS = [
    (
        "stockholm_case55_clean/building03_domain.pslg",
        (273, 188, 10),
        "invalid PSLG input",
    ),
    (
        "stockholm_case55_clean/building17_domain.pslg",
        (156, 103, 4),
        "invalid PSLG input",
    ),
    (
        "stockholm_case55_clean/building20_domain.pslg",
        (366, 233, 11),
        "invalid PSLG input",
    ),
    (
        "stockholm_case55_clean/building31_domain.pslg",
        (84, 61, 0),
        "invalid mesh topology",
    ),
]


def _case_path(name: str) -> Path:
    case_dir = Path(__file__).resolve().parents[1] / "cases"
    return case_dir / name


@pytest.mark.parametrize(
    ("name", "expected_sizes", "_expected_error"),
    _CASE55_BAD_DOMAINS,
)
def test_stockholm_case55_domain_files_round_trip_shapes(
    name: str,
    expected_sizes: tuple[int, int, int],
    _expected_error: str,
):
    points, segments, holes = dm.read_domain_file(_case_path(name))

    assert points.shape == (expected_sizes[0], 2)
    assert segments is not None
    assert segments.shape == (expected_sizes[1], 2)
    if expected_sizes[2] == 0:
        assert holes is None
    else:
        assert holes is not None
        assert holes.shape == (expected_sizes[2], 2)


@pytest.mark.parametrize(
    ("name", "_expected_sizes", "expected_error"),
    _CASE55_BAD_DOMAINS,
)
def test_stockholm_case55_domains_capture_current_failure_modes(
    name: str,
    _expected_sizes: tuple[int, int, int],
    expected_error: str,
):
    with pytest.raises(RuntimeError, match=re.escape(expected_error)):
        dm.generate_file(_case_path(name), min_angle=25.0, refine=False)
