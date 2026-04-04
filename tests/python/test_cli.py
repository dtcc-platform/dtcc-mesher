from __future__ import annotations

import math
import os
import subprocess
from pathlib import Path

import pytest

import dtcc_mesher as dm


def _parse_summary_text(text: str) -> dict[str, float]:
    summary: dict[str, float] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        summary[key] = float(value)
    return summary


def test_cli_and_python_match(tmp_path):
    cli = os.environ.get("DTCC_MESHER_TEST_CLI")
    c_api_dump = os.environ.get("DTCC_MESHER_TEST_C_API_DUMP")
    case_dir = os.environ.get("DTCC_MESHER_TEST_CASE_DIR")

    if not cli or not c_api_dump or not case_dir:
        pytest.skip("CLI consistency paths are not configured")

    input_path = Path(case_dir) / "square_hole_domain.pslg"
    out_base = tmp_path / "square_hole_cli"

    subprocess.run([cli, str(input_path), str(out_base)], check=True)
    cli_summary = _parse_summary_text((tmp_path / "square_hole_cli.summary.txt").read_text(encoding="utf-8"))

    mesh = dm.mesh(dm.read_domain(input_path))
    c_api_summary = _parse_summary_text(
        subprocess.run([c_api_dump, str(input_path)], check=True, capture_output=True, text=True).stdout
    )

    assert int(cli_summary["triangle_count"]) == mesh.summary.triangle_count == int(c_api_summary["triangle_count"])
    assert int(cli_summary["point_count"]) == mesh.summary.point_count == int(c_api_summary["point_count"])
    assert int(cli_summary["count_min_angle_lt_20"]) == mesh.summary.count_min_angle_lt_20 == int(
        c_api_summary["count_min_angle_lt_20"]
    )
    assert math.isclose(cli_summary["area_min"], mesh.summary.area_min, rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(cli_summary["area_max"], mesh.summary.area_max, rel_tol=1e-9, abs_tol=1e-9)
