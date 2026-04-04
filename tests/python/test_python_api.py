from __future__ import annotations

import os
import subprocess
import sys
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


def test_generate_pslg_without_holes():
    mesh = dm.generate_file(_case_path("square_domain.pslg"))
    assert mesh.triangles.shape[1] == 3
    assert mesh.summary.triangle_count > 0


def test_generate_pslg_max_edge_length_increases_density_without_refinement():
    coarse = dm.generate_file(_case_path("square_domain.pslg"), refine=False)
    split = dm.generate_file(_case_path("square_domain.pslg"), refine=False, max_edge_length=0.5)

    assert coarse.summary.triangle_count == 2
    assert split.summary.triangle_count > coarse.summary.triangle_count


def test_generate_pslg_max_area_refines_even_when_refine_flag_is_disabled():
    coarse = dm.generate_file(_case_path("square_domain.pslg"), refine=False)
    refined = dm.generate_file(_case_path("square_domain.pslg"), refine=False, max_area=0.1)

    assert coarse.summary.triangle_count == 2
    assert refined.summary.triangle_count > coarse.summary.triangle_count
    assert refined.summary.area_max <= 0.1 + 1e-8


def test_generate_rejects_negative_max_area():
    with pytest.raises(RuntimeError, match="max_area must be non-negative"):
        dm.generate_file(_case_path("square_domain.pslg"), max_area=-1.0)


def test_generate_graph_returns_region_markers_for_shared_edge_partition():
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [0.5, 1.0],
        ],
        dtype=np.float64,
    )
    segments = np.array(
        [
            [0, 4],
            [4, 1],
            [1, 2],
            [2, 5],
            [5, 3],
            [3, 0],
            [4, 5],
        ],
        dtype=np.uint32,
    )
    region_points = np.array([[0.25, 0.5], [0.75, 0.5]], dtype=np.float64)
    region_markers = np.array([10, 20], dtype=np.int32)

    mesh = dm.generate_graph(
        points,
        segments=segments,
        region_points=region_points,
        region_markers=region_markers,
        max_edge_length=0.25,
        refine=False,
    )

    assert mesh.markers is not None
    assert set(np.asarray(mesh.markers, dtype=int)) == {10, 20}

    centroids = mesh.points[mesh.triangles].mean(axis=1)
    for marker, centroid in zip(np.asarray(mesh.markers, dtype=int), centroids):
        if centroid[0] < 0.5:
            assert marker == 10
        else:
            assert marker == 20


def test_generate_coverage_builds_single_region_graph_from_polygons():
    shapely = pytest.importorskip("shapely.geometry")
    left = shapely.box(0.0, 0.0, 0.5, 1.0)
    right = shapely.box(0.5, 0.0, 1.0, 1.0)

    mesh = dm.generate_coverage(
        [left, right],
        markers=[10, 20],
        max_edge_length=0.25,
        refine=False,
    )

    assert mesh.markers is not None
    assert set(np.asarray(mesh.markers, dtype=int)) == {10, 20}
    assert mesh.summary.triangle_count > 2


def test_generate_coverage_preserves_closed_outer_boundary_sampling():
    shapely = pytest.importorskip("shapely.geometry")
    mesh = dm.generate_coverage(
        [shapely.box(0.0, 0.0, 500.0, 500.0)],
        markers=[1],
        max_edge_length=10.0,
        refine=False,
    )

    bottom_vertices = np.count_nonzero(np.abs(mesh.points[:, 1]) < 1e-9)
    assert bottom_vertices > 10

    aspect_ratios = []
    for triangle in mesh.triangles:
        points = mesh.points[triangle]
        edge_lengths = np.array(
            [
                np.linalg.norm(points[1] - points[0]),
                np.linalg.norm(points[2] - points[1]),
                np.linalg.norm(points[0] - points[2]),
            ]
        )
        longest = float(edge_lengths.max())
        area = float(abs(np.cross(points[1] - points[0], points[2] - points[0])) / 2.0)
        altitude = (2.0 * area / longest) if longest > 0.0 else 0.0
        aspect_ratios.append(longest / altitude if altitude > 0.0 else np.inf)

    assert max(aspect_ratios) < 20.0


def test_generate_graph_reports_t_junctions_in_coverage_input():
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [0.5, 0.5],
        ],
        dtype=np.float64,
    )
    segments = np.array(
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [4, 5],
        ],
        dtype=np.uint32,
    )
    region_points = np.array([[0.25, 0.25]], dtype=np.float64)
    region_markers = np.array([1], dtype=np.int32)

    with pytest.raises(RuntimeError, match="endpoint lies on another segment interior"):
        dm.generate_graph(
            points,
            segments=segments,
            region_points=region_points,
            region_markers=region_markers,
            refine=False,
        )


def test_plot_mesh_with_summary_smoke():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mesh = dm.generate(
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, 0.5),
        ]
    )

    fig, (ax_mesh, ax_summary) = dm.plot_mesh_with_summary(mesh, title="square_center")
    assert ax_mesh.get_title() == "square_center"
    assert len(ax_summary.texts) == 1
    assert "triangles:" in ax_summary.texts[0].get_text()
    plt.close(fig)


def test_invalid_points_shape():
    with pytest.raises(ValueError):
        dm.generate([(0.0, 0.0, 1.0)])


def test_plot_demo_smoke():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    repo_root = Path(__file__).resolve().parents[2]
    demo_path = repo_root / "demos" / "python" / "plot_mesh.py"

    result = subprocess.run(
        [sys.executable, str(demo_path), "square_hole_domain", "--quiet"],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )

    assert "summary:" in result.stdout
    assert "triangles:" in result.stdout
