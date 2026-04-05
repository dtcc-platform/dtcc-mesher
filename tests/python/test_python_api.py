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


def _mesh_domain_file(path: Path, **option_kwargs) -> dm.Mesh:
    return dm.mesh(
        dm.read_domain(path),
        options=dm.MeshingOptions(**option_kwargs),
    )


def test_mesh_points():
    mesh = dm.mesh(
        dm.Domain(
            points=[
                (0.0, 0.0),
                (1.0, 0.0),
                (1.0, 1.0),
                (0.0, 1.0),
                (0.5, 0.5),
            ]
        )
    )

    assert mesh.points.shape == (5, 2)
    assert mesh.triangles.shape == (4, 3)
    assert mesh.points.dtype == np.float64
    assert mesh.triangles.dtype == np.uint32
    assert mesh.summary.triangle_count == 4
    assert mesh.summary.count_min_angle_lt_20 == 0


def test_mesh_domain_object():
    domain = dm.Domain(
        points=[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, 0.5),
        ]
    )

    mesh = dm.mesh(domain)

    assert mesh.points.shape == (5, 2)
    assert mesh.triangles.shape == (4, 3)
    assert mesh.summary.triangle_count == 4


def test_read_domain():
    domain = dm.read_domain(_case_path("square_domain.pslg"))

    assert domain.points.shape == (4, 2)
    assert domain.segments is not None
    assert domain.segments.shape == (4, 2)
    assert domain.holes is not None
    assert domain.holes.shape == (0, 2)


def test_domain_from_loops_with_hole():
    domain = dm.Domain.from_loops(
        [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
        holes=[[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]],
    )

    mesh = dm.mesh(domain, options=dm.MeshingOptions(refine=False))

    assert domain.segments is not None
    assert domain.holes is not None
    assert domain.holes.shape == (1, 2)
    assert mesh.summary.triangle_count > 0


def test_mesh_domain_file_and_write_svg(tmp_path):
    mesh = _mesh_domain_file(_case_path("square_hole_domain.pslg"))
    svg_path = tmp_path / "mesh.svg"
    mesh.write_svg(svg_path)
    assert svg_path.exists()
    assert svg_path.read_text(encoding="utf-8").startswith("<svg")


def test_mesh_domain_without_holes():
    mesh = _mesh_domain_file(_case_path("square_domain.pslg"))
    assert mesh.triangles.shape[1] == 3
    assert mesh.summary.triangle_count > 0


def test_mesh_domain_max_edge_length_increases_density_without_refinement():
    coarse = _mesh_domain_file(_case_path("square_domain.pslg"), refine=False)
    split = _mesh_domain_file(_case_path("square_domain.pslg"), refine=False, max_edge_length=0.5)

    assert coarse.summary.triangle_count == 2
    assert split.summary.triangle_count > coarse.summary.triangle_count


def test_mesh_domain_max_area_refines_even_when_refine_flag_is_disabled():
    coarse = _mesh_domain_file(_case_path("square_domain.pslg"), refine=False)
    refined = _mesh_domain_file(_case_path("square_domain.pslg"), refine=False, max_area=0.1)

    assert coarse.summary.triangle_count == 2
    assert refined.summary.triangle_count > coarse.summary.triangle_count
    assert refined.summary.area_max <= 0.1 + 1e-8


def test_mesh_domain_rejects_negative_max_area():
    with pytest.raises(RuntimeError, match="max_area must be non-negative"):
        _mesh_domain_file(_case_path("square_domain.pslg"), max_area=-1.0)


def test_mesh_graph_returns_region_markers_for_shared_edge_partition():
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

    mesh = dm.mesh(
        dm.CoverageGraph(
            points=points,
            segments=segments,
            region_points=region_points,
            region_markers=region_markers,
        ),
        options=dm.MeshingOptions(max_edge_length=0.25, refine=False),
    )

    assert mesh.markers is not None
    assert set(np.asarray(mesh.markers, dtype=int)) == {10, 20}

    centroids = mesh.points[mesh.triangles].mean(axis=1)
    for marker, centroid in zip(np.asarray(mesh.markers, dtype=int), centroids):
        if centroid[0] < 0.5:
            assert marker == 10
        else:
            assert marker == 20


def test_mesh_coverage_object_builds_single_region_graph_from_polygons():
    shapely = pytest.importorskip("shapely.geometry")
    left = shapely.box(0.0, 0.0, 0.5, 1.0)
    right = shapely.box(0.5, 0.0, 1.0, 1.0)

    mesh = dm.mesh(
        dm.Coverage([left, right], markers=[10, 20]),
        options=dm.MeshingOptions(max_edge_length=0.25, refine=False),
    )

    assert mesh.markers is not None
    assert set(np.asarray(mesh.markers, dtype=int)) == {10, 20}
    assert mesh.summary.triangle_count > 2


def test_mesh_coverage_preserves_closed_outer_boundary_sampling():
    shapely = pytest.importorskip("shapely.geometry")
    mesh = dm.mesh(
        dm.Coverage([shapely.box(0.0, 0.0, 500.0, 500.0)], markers=[1]),
        options=dm.MeshingOptions(max_edge_length=10.0, refine=False),
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


def test_mesh_coverage_matches_domain_for_large_coordinate_acute_polygon():
    shapely = pytest.importorskip("shapely.geometry")
    polygon = shapely.Polygon(
        [
            (675000.0, 6581000.0),
            (675020.0, 6581000.0),
            (675001.0, 6581004.0),
            (675000.0, 6581020.0),
            (675000.0, 6581000.0),
        ]
    )

    coverage_mesh = dm.mesh(
        dm.Coverage([polygon], markers=[1]),
        options=dm.MeshingOptions(min_angle=25.0, refine=True, max_refine_steps=5000),
    )
    domain_mesh = dm.mesh(
        dm.Domain.from_loops(polygon.exterior.coords),
        options=dm.MeshingOptions(min_angle=25.0, refine=True, max_refine_steps=5000),
    )

    assert coverage_mesh.triangles.shape[0] == domain_mesh.triangles.shape[0]
    assert coverage_mesh.markers is not None
    assert set(np.asarray(coverage_mesh.markers, dtype=int)) == {1}


def test_mesh_coverage_default_tolerance_merges_city_scale_micro_edges():
    shapely = pytest.importorskip("shapely.geometry")
    polygon = shapely.Polygon(
        [
            (675000.0, 6580000.0),
            (675100.0, 6580000.0),
            (675100.0, 6580100.0),
            (675000.0001, 6580100.00005),
            (675000.0, 6580100.0),
        ]
    )

    default_graph = dm.Coverage([polygon], [1]).graph(max_edge_length=10.0)
    exact_graph = dm.Coverage([polygon], [1], tolerance=1e-9).graph(max_edge_length=10.0)

    default_lengths = np.linalg.norm(
        default_graph.points[default_graph.segments[:, 1]] - default_graph.points[default_graph.segments[:, 0]],
        axis=1,
    )
    exact_lengths = np.linalg.norm(
        exact_graph.points[exact_graph.segments[:, 1]] - exact_graph.points[exact_graph.segments[:, 0]],
        axis=1,
    )

    assert default_graph.segments.shape[0] == 4
    assert exact_graph.segments.shape[0] == 5
    assert float(default_lengths.min()) > 1.0
    assert float(exact_lengths.min()) < 1e-3


def test_mesh_coverage_handles_courtyard_region_partition():
    shapely = pytest.importorskip("shapely.geometry")

    ground = shapely.Polygon(
        [(0.0, 0.0), (12.0, 0.0), (12.0, 12.0), (0.0, 12.0), (0.0, 0.0)],
        holes=[[(3.0, 3.0), (9.0, 3.0), (9.0, 9.0), (3.0, 9.0), (3.0, 3.0)]],
    )
    building = shapely.Polygon(
        [(3.0, 3.0), (9.0, 3.0), (9.0, 9.0), (3.0, 9.0), (3.0, 3.0)],
        holes=[[(4.5, 4.5), (7.5, 4.5), (7.5, 7.5), (4.5, 7.5), (4.5, 4.5)]],
    )
    courtyard = shapely.box(4.5, 4.5, 7.5, 7.5)

    mesh = dm.mesh(
        dm.Coverage([ground, building, courtyard], markers=[-2, 10, -2]),
        options=dm.MeshingOptions(max_edge_length=1.5, min_angle=25.0, refine=False),
    )

    assert mesh.markers is not None
    assert set(np.asarray(mesh.markers, dtype=int)) == {-2, 10}
    assert mesh.summary.triangle_count > 0


def test_mesh_graph_reports_t_junctions_in_coverage_input():
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
        dm.mesh(
            dm.CoverageGraph(
                points=points,
                segments=segments,
                region_points=region_points,
                region_markers=region_markers,
            ),
            options=dm.MeshingOptions(refine=False),
        )


def test_mesh_graph_handles_stockholm_case15_recovery_regression():
    with np.load(_case_path("stockholm_case15_graph_maxh10.npz")) as data:
        graph = dm.CoverageGraph(
            points=data["points"],
            segments=data["segments"],
            region_points=data["region_points"],
            region_markers=data["region_markers"],
        )

    mesh = dm.mesh(
        graph,
        options=dm.MeshingOptions(min_angle=25.0, max_edge_length=10.0, refine=False),
    )

    assert mesh.summary.triangle_count > 7000
    assert mesh.markers is not None
    assert len(mesh.markers) == mesh.summary.triangle_count


def test_mesh_defaults_protect_33_degree_input_corner():
    length = 10.0
    angle = np.deg2rad(33.0)
    mesh = dm.mesh(
        dm.Domain.from_loops(
            [
                (0.0, 0.0),
                (length, 0.0),
                (length * np.cos(angle), length * np.sin(angle)),
                (0.0, 0.0),
            ]
        ),
        options=dm.MeshingOptions(min_angle=25.0, refine=True),
    )

    assert mesh.summary.protected_corner_count == 1
    assert mesh.summary.triangle_count > 1


def test_mesh_graph_case55_unrestricted_reports_step_limit_cleanly():
    with np.load(_case_path("stockholm_case55_graph_unrestricted.npz")) as data:
        graph = dm.CoverageGraph(
            points=data["points"],
            segments=data["segments"],
            region_points=data["region_points"],
            region_markers=data["region_markers"],
        )

    with pytest.raises(RuntimeError, match=r"quality refinement reached step limit after 50 steps"):
        dm.mesh(
            graph,
            options=dm.MeshingOptions(min_angle=25.0, refine=True, max_refine_steps=50),
        )


def test_mesh_graph_case55_unrestricted_completes_with_good_quality():
    with np.load(_case_path("stockholm_case55_graph_unrestricted.npz")) as data:
        graph = dm.CoverageGraph(
            points=data["points"],
            segments=data["segments"],
            region_points=data["region_points"],
            region_markers=data["region_markers"],
        )

    mesh = dm.mesh(
        graph,
        options=dm.MeshingOptions(min_angle=25.0, refine=True, max_refine_steps=2000),
    )

    assert mesh.summary.triangle_count > 4500
    assert mesh.summary.radius_edge_ratio_max < 1.6
    assert mesh.summary.count_min_angle_lt_20 <= 2


def test_mesh_graph_case34_unrestricted_completes_with_good_quality():
    with np.load(_case_path("stockholm_case34_graph_unrestricted.npz")) as data:
        graph = dm.CoverageGraph(
            points=data["points"],
            segments=data["segments"],
            region_points=data["region_points"],
            region_markers=data["region_markers"],
        )

    mesh = dm.mesh(
        graph,
        options=dm.MeshingOptions(min_angle=25.0, refine=True, max_refine_steps=20000),
    )

    assert mesh.summary.triangle_count > 7500
    assert mesh.summary.radius_edge_ratio_max < 1.6
    assert mesh.summary.count_min_angle_lt_20 <= 2


def test_plot_mesh_with_summary_smoke():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mesh = dm.mesh(
        dm.Domain(
            points=[
                (0.0, 0.0),
                (1.0, 0.0),
                (1.0, 1.0),
                (0.0, 1.0),
                (0.5, 0.5),
            ]
        )
    )

    fig, (ax_mesh, ax_summary) = dm.plot_mesh_with_summary(mesh, title="square_center")
    assert ax_mesh.get_title() == "square_center"
    assert len(ax_summary.texts) == 1
    assert "triangles:" in ax_summary.texts[0].get_text()
    plt.close(fig)


def test_invalid_points_shape():
    with pytest.raises(ValueError):
        dm.mesh(dm.Domain(points=[(0.0, 0.0, 1.0)]))


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
