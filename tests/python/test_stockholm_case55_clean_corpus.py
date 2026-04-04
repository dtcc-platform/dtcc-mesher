from __future__ import annotations

import json
from pathlib import Path

import dtcc_mesher as dm


def _case55_clean_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "cases" / "stockholm_case55_clean"


def _manifest() -> dict:
    return json.loads((_case55_clean_dir() / "manifest.json").read_text(encoding="utf-8"))


def test_stockholm_case55_clean_manifest_metadata():
    manifest = _manifest()
    summary = manifest["summary"]

    assert manifest["source_case"] == 55
    assert manifest["crs"] == "EPSG:3006"
    assert manifest["mesher_input"] == {
        "max_mesh_size": 10.0,
        "min_angle": 25.0,
        "refine": False,
    }
    assert summary["domain_count"] == 89
    assert summary["ground_domain_count"] == 55
    assert summary["building_domain_count"] == 34
    assert summary["failing_domain_count"] == 3
    assert len(manifest["domains"]) == summary["domain_count"]


def test_stockholm_case55_clean_domain_files_match_manifest():
    case_dir = _case55_clean_dir()
    manifest = _manifest()

    for entry in manifest["domains"]:
        path = case_dir / entry["filename"]
        assert path.exists()

        points, segments, holes = dm.read_domain_file(path)

        assert points.shape == (entry["point_count"], 2)
        assert segments is not None
        assert segments.shape == (entry["segment_count"], 2)
        if entry["hole_marker_count"] == 0:
            assert holes is None
        else:
            assert holes is not None
            assert holes.shape == (entry["hole_marker_count"], 2)


def test_stockholm_case55_clean_generate_partition_matches_manifest():
    case_dir = _case55_clean_dir()
    manifest = _manifest()
    min_angle = manifest["mesher_input"]["min_angle"]
    refine = manifest["mesher_input"]["refine"]

    seen_failures: dict[str, str] = {}
    success_count = 0

    for entry in manifest["domains"]:
        path = case_dir / entry["filename"]

        try:
            dm.generate_file(path, min_angle=min_angle, refine=refine)
        except RuntimeError as exc:
            seen_failures[entry["filename"]] = str(exc)
            assert entry["status"] == "error"
            assert entry["error"] in str(exc)
        else:
            success_count += 1
            assert entry["status"] == "ok"

    expected_failures = {
        entry["filename"]: entry["error"]
        for entry in manifest["domains"]
        if entry["status"] == "error"
    }

    assert success_count == manifest["summary"]["domain_count"] - manifest["summary"]["failing_domain_count"]
    assert seen_failures == expected_failures
