from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import dtcc_mesher as dm


def _case_path(name: str) -> Path:
    case_dir = os.environ.get("DTCC_MESHER_TEST_CASE_DIR")
    if case_dir:
        return Path(case_dir) / name
    return Path(__file__).resolve().parents[1] / "cases" / name


def test_generate_points():
    mesh = dm.generate(
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, 0.5),
        ]
    )

    assert mesh.points.shape == (5, 2)
    assert mesh.triangles.shape == (4, 3)
    assert mesh.points.dtype == np.float64
    assert mesh.triangles.dtype == np.uint32
    assert mesh.summary.triangle_count == 4
    assert mesh.summary.count_min_angle_lt_20 == 0


def test_generate_pslg_and_write_svg(tmp_path):
    mesh = dm.generate_file(_case_path("square_hole_domain.pslg"))
    svg_path = tmp_path / "mesh.svg"
    mesh.write_svg(svg_path)
    assert svg_path.exists()
    assert svg_path.read_text(encoding="utf-8").startswith("<svg")


def test_invalid_points_shape():
    with pytest.raises(ValueError):
        dm.generate([(0.0, 0.0, 1.0)])
