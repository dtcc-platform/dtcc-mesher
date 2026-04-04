from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from . import _core
from . import _io
from ._coverage import build_coverage_domain
from ._plot import format_summary_lines, plot_mesh, plot_mesh_with_summary, show_mesh

__all__ = [
    "Mesh",
    "QualitySummary",
    "generate",
    "generate_coverage",
    "generate_graph",
    "generate_file",
    "format_summary_lines",
    "plot_mesh",
    "plot_mesh_with_summary",
    "read_domain_file",
    "show_mesh",
]

__version__ = "0.1.0"


@dataclass(slots=True)
class QualitySummary:
    point_count: int
    input_point_count: int
    steiner_point_count: int
    segment_split_point_count: int
    triangle_split_point_count: int
    protected_corner_count: int
    exempt_triangle_count: int
    triangle_count: int
    area_min: float
    area_mean: float
    area_max: float
    min_angle_deg_min: float
    min_angle_deg_mean: float
    min_angle_deg_max: float
    edge_ratio_min: float
    edge_ratio_mean: float
    edge_ratio_max: float
    radius_edge_ratio_min: float
    radius_edge_ratio_mean: float
    radius_edge_ratio_max: float
    count_min_angle_lt_20: int
    count_min_angle_lt_30: int


@dataclass
class Mesh:
    points: np.ndarray
    triangles: np.ndarray
    segments: np.ndarray
    summary: QualitySummary
    markers: np.ndarray | None = None

    def write_svg(self, path: str | Path) -> None:
        _io.write_svg(self, path)

    def write_triangles(self, path: str | Path) -> None:
        _io.write_triangles(self, path)

    def write_quality_csv(self, path: str | Path) -> None:
        _io.write_quality_csv(self, path)

    def write_quality_summary(self, path: str | Path) -> None:
        _io.write_quality_summary(self, path)

    def plot(self, **kwargs):
        return plot_mesh(self, **kwargs)

    def show(self, **kwargs):
        return show_mesh(self, **kwargs)

    def summary_lines(self) -> list[str]:
        return format_summary_lines(self.summary)

    def analyze(self) -> QualitySummary:
        raw = _core._analyze_raw(
            self.points,
            self.triangles,
            self.segments,
            self.summary.input_point_count,
            self.summary.segment_split_point_count,
            self.summary.triangle_split_point_count,
            self.summary.protected_corner_count,
            self.summary.exempt_triangle_count,
        )
        self.summary = QualitySummary(**raw)
        return self.summary


def _coerce_points(points: object | None) -> np.ndarray | None:
    if points is None:
        return None
    array = np.asarray(points, dtype=np.float64)
    if array.size == 0:
        return np.empty((0, 2), dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("expected an array with shape (N, 2)")
    return np.ascontiguousarray(array)


def _coerce_segments(segments: object | None) -> np.ndarray | None:
    if segments is None:
        return None
    array = np.asarray(segments, dtype=np.uint32)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("expected an array with shape (M, 2)")
    return np.ascontiguousarray(array)


def _acute_mode_value(mode: Literal["none", "simple", "shell"]) -> int:
    if mode == "none":
        return 0
    if mode == "simple":
        return 1
    if mode == "shell":
        return 2
    raise ValueError("acute_protection must be 'none', 'simple', or 'shell'")


def _summary_from_raw(raw: dict) -> QualitySummary:
    return QualitySummary(**raw)


def _coerce_markers(markers: object | None) -> np.ndarray | None:
    if markers is None:
        return None
    array = np.asarray(markers, dtype=np.int32)
    if array.ndim != 1:
        raise ValueError("expected a marker array with shape (N,)")
    return np.ascontiguousarray(array)


def generate(
    points: object,
    *,
    segments: object | None = None,
    holes: object | None = None,
    min_angle: float = 20.0,
    max_area: float | None = None,
    max_edge_length: float | None = None,
    refine: bool = True,
    off_centers: bool = False,
    verbose: bool = False,
    acute_protection: Literal["none", "simple", "shell"] = "shell",
    protect_angle: float | None = None,
    max_refine_steps: int = 0,
    max_protection_levels: int = 6,
) -> Mesh:
    points_array = _coerce_points(points)
    segments_array = _coerce_segments(segments)
    holes_array = _coerce_points(holes)

    raw = _core._generate_raw(
        points_array,
        segments_array,
        holes_array,
        min_angle,
        max_area,
        max_edge_length,
        refine,
        off_centers,
        verbose,
        _acute_mode_value(acute_protection),
        protect_angle,
        max_refine_steps,
        max_protection_levels,
    )

    return Mesh(
        points=np.asarray(raw["points"], dtype=np.float64),
        triangles=np.asarray(raw["triangles"], dtype=np.uint32),
        segments=np.asarray(raw["segments"], dtype=np.uint32),
        summary=_summary_from_raw(raw["summary"]),
        markers=np.asarray(raw["markers"], dtype=np.int32) if "markers" in raw else None,
    )


def generate_graph(
    points: object,
    *,
    segments: object,
    region_points: object,
    region_markers: object,
    min_angle: float = 20.0,
    max_area: float | None = None,
    max_edge_length: float | None = None,
    refine: bool = True,
    off_centers: bool = False,
    verbose: bool = False,
    acute_protection: Literal["none", "simple", "shell"] = "shell",
    protect_angle: float | None = None,
    max_refine_steps: int = 0,
    max_protection_levels: int = 6,
) -> Mesh:
    points_array = _coerce_points(points)
    segments_array = _coerce_segments(segments)
    region_points_array = _coerce_points(region_points)
    region_markers_array = _coerce_markers(region_markers)

    if points_array is None or segments_array is None:
        raise ValueError("points and segments are required for coverage meshing")
    if region_points_array is None or region_markers_array is None:
        raise ValueError("region_points and region_markers are required for coverage meshing")

    raw = _core._generate_coverage_raw(
        points_array,
        segments_array,
        region_points_array,
        region_markers_array,
        min_angle,
        max_area,
        max_edge_length,
        refine,
        off_centers,
        verbose,
        _acute_mode_value(acute_protection),
        protect_angle,
        max_refine_steps,
        max_protection_levels,
    )

    return Mesh(
        points=np.asarray(raw["points"], dtype=np.float64),
        triangles=np.asarray(raw["triangles"], dtype=np.uint32),
        segments=np.asarray(raw["segments"], dtype=np.uint32),
        summary=_summary_from_raw(raw["summary"]),
        markers=np.asarray(raw["markers"], dtype=np.int32),
    )


def generate_coverage(
    polygons: object,
    *,
    markers: object,
    min_angle: float = 20.0,
    max_area: float | None = None,
    max_edge_length: float | None = None,
    refine: bool = True,
    off_centers: bool = False,
    verbose: bool = False,
    acute_protection: Literal["none", "simple", "shell"] = "shell",
    protect_angle: float | None = None,
    max_refine_steps: int = 0,
    max_protection_levels: int = 6,
    tolerance: float = 1e-9,
) -> Mesh:
    coverage_points, coverage_segments, region_points, region_markers = build_coverage_domain(
        polygons,
        markers,
        max_edge_length=max_edge_length,
        tolerance=tolerance,
    )

    return generate_graph(
        coverage_points,
        segments=coverage_segments,
        region_points=region_points,
        region_markers=region_markers,
        min_angle=min_angle,
        max_area=max_area,
        max_edge_length=max_edge_length,
        refine=refine,
        off_centers=off_centers,
        verbose=verbose,
        acute_protection=acute_protection,
        protect_angle=protect_angle,
        max_refine_steps=max_refine_steps,
        max_protection_levels=max_protection_levels,
    )


def read_domain_file(path: str | Path) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    return _io.read_domain_file(path)


def generate_file(path: str | Path, **kwargs) -> Mesh:
    points, segments, holes = read_domain_file(path)
    return generate(points, segments=segments, holes=holes, **kwargs)
