from __future__ import annotations

import re
from pathlib import Path

import pytest

import dtcc_mesher as dm


_CASE55_BAD_DOMAINS = [
    (
        "stockholm_case55_clean/building03_domain.pslg",
        (273, 188, 10),
        "segment graph endpoint lies on another segment interior",
    ),
    (
        "stockholm_case55_clean/building17_domain.pslg",
        (156, 103, 4),
        "segment graph endpoint lies on another segment interior",
    ),
    (
        "stockholm_case55_clean/building20_domain.pslg",
        (366, 233, 11),
        "segment graph endpoint lies on another segment interior",
    ),
]

_CASE55_FIXED_DOMAIN = "stockholm_case55_clean/building31_domain.pslg"


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
    domain = dm.read_domain(_case_path(name))
    points = domain.points
    segments = domain.segments
    holes = domain.holes

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
        dm.mesh(
            dm.read_domain(_case_path(name)),
            options=dm.MeshingOptions(min_angle=25.0, refine=False),
        )


def test_stockholm_case55_building31_domain_now_succeeds():
    mesh = dm.mesh(
        dm.read_domain(_case_path(_CASE55_FIXED_DOMAIN)),
        options=dm.MeshingOptions(min_angle=25.0, refine=False),
    )

    assert mesh.points.shape == (85, 2)
    assert mesh.triangles.shape == (106, 3)
    assert mesh.segments.shape == (62, 2)
