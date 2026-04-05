from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Literal

import numpy as np

from . import _core
from . import _io
from ._coverage import build_coverage_domain
from ._plot import format_summary_lines, plot_mesh, plot_mesh_with_summary, show_mesh

__all__ = [
    "Coverage",
    "CoverageGraph",
    "Domain",
    "Mesh",
    "MeshingOptions",
    "QualitySummary",
    "format_summary_lines",
    "mesh",
    "plot_mesh",
    "plot_mesh_with_summary",
    "read_domain",
    "show_mesh",
]

__version__ = "0.1.0"


@dataclass(slots=True)
class MeshingOptions:
    """Meshing controls for :func:`mesh`.

    ``max_edge_length`` is the preferred public size control for 2D city
    meshing. It represents the maximum target constrained-edge length and
    nominal mesh spacing in the input plane. Use ``None`` to leave the mesh
    globally unrestricted and let geometry plus ``min_angle`` drive
    refinement.

    ``max_area`` is a lower-level area cap for backends that refine by area.
    Most callers should prefer ``max_edge_length``.

    Refinement uses a standard constrained-Delaunay scheme: encroached
    constrained segments are split at midpoints and bad triangles are split
    at circumcenters.

    With ``acute_protection="shell"`` (the default), PSLG corners below
    60 degrees are protected before refinement unless ``protect_angle`` is
    set explicitly.
    """

    min_angle: float = 20.0
    max_area: float | None = None
    max_edge_length: float | None = None
    refine: bool = True
    off_centers: bool = False
    verbose: bool = False
    acute_protection: Literal["none", "simple", "shell"] = "shell"
    protect_angle: float | None = None
    max_refine_steps: int = 0
    max_protection_levels: int = 6


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


def _require_points(points: object, *, name: str) -> np.ndarray:
    array = _coerce_points(points)
    if array is None:
        raise ValueError(f"{name} are required")
    return array


def _coerce_segments(segments: object | None) -> np.ndarray | None:
    if segments is None:
        return None
    array = np.asarray(segments, dtype=np.uint32)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("expected an array with shape (M, 2)")
    return np.ascontiguousarray(array)


def _require_segments(segments: object, *, name: str) -> np.ndarray:
    array = _coerce_segments(segments)
    if array is None:
        raise ValueError(f"{name} are required")
    return array


def _coerce_markers(markers: object | None) -> np.ndarray | None:
    if markers is None:
        return None
    array = np.asarray(markers, dtype=np.int32)
    if array.ndim != 1:
        raise ValueError("expected a marker array with shape (N,)")
    return np.ascontiguousarray(array)


def _require_markers(markers: object, *, name: str) -> np.ndarray:
    array = _coerce_markers(markers)
    if array is None:
        raise ValueError(f"{name} are required")
    return array


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


def _mesh_from_raw(raw: dict) -> Mesh:
    return Mesh(
        points=np.asarray(raw["points"], dtype=np.float64),
        triangles=np.asarray(raw["triangles"], dtype=np.uint32),
        segments=np.asarray(raw["segments"], dtype=np.uint32),
        summary=_summary_from_raw(raw["summary"]),
        markers=np.asarray(raw["markers"], dtype=np.int32) if "markers" in raw else None,
    )


def _sanitize_loop(points: object, tol: float = 1e-10) -> np.ndarray:
    array = _require_points(points, name="loop")
    clean: list[np.ndarray] = []
    for point in array:
        if not clean or np.linalg.norm(point - clean[-1]) > tol:
            clean.append(point)

    if len(clean) > 1 and np.linalg.norm(clean[0] - clean[-1]) <= tol:
        clean.pop()

    clean_array = np.asarray(clean, dtype=np.float64)
    if len(clean_array) < 3:
        raise ValueError("expected a loop with at least three distinct vertices")
    return clean_array


def _point_in_ring(point: np.ndarray, ring: np.ndarray) -> bool:
    x = float(point[0])
    y = float(point[1])
    inside = False

    for start, end in zip(ring, np.vstack([ring[1:], ring[:1]])):
        x0 = float(start[0])
        y0 = float(start[1])
        x1 = float(end[0])
        y1 = float(end[1])
        intersects = (y0 > y) != (y1 > y)
        if not intersects:
            continue
        cross_x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
        if x < cross_x:
            inside = not inside

    return inside


def _loop_interior_point(loop: np.ndarray) -> np.ndarray:
    min_y = float(np.min(loop[:, 1]))
    max_y = float(np.max(loop[:, 1]))
    span_y = max_y - min_y
    if span_y <= 0.0:
        raise ValueError("loop must span a positive area")

    sample_offsets = (0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875)
    edges = np.vstack([loop[1:], loop[:1]])

    for offset in sample_offsets:
        y = min_y + offset * span_y
        intersections: list[float] = []
        for start, end in zip(loop, edges):
            y0 = float(start[1])
            y1 = float(end[1])
            if abs(y1 - y0) <= 1e-12:
                continue
            low_y = min(y0, y1)
            high_y = max(y0, y1)
            if not (low_y <= y < high_y):
                continue
            t = (y - y0) / (y1 - y0)
            x = float(start[0]) + t * (float(end[0]) - float(start[0]))
            intersections.append(x)

        intersections.sort()
        for left, right in zip(intersections[0::2], intersections[1::2]):
            if right - left <= 1e-12:
                continue
            candidate = np.array([(left + right) * 0.5, y], dtype=np.float64)
            if _point_in_ring(candidate, loop):
                return candidate

    centroid = np.mean(loop, axis=0, dtype=np.float64)
    if _point_in_ring(centroid, loop):
        return centroid

    raise ValueError("could not compute an interior point for loop")


class _DomainBuilder:
    def __init__(self, tolerance: float = 1e-10):
        self._scale = 1.0 / tolerance
        self._points: list[np.ndarray] = []
        self._segments: list[tuple[int, int]] = []
        self._point_index: dict[tuple[int, int], int] = {}
        self._segment_index: set[tuple[int, int]] = set()

    def _key(self, point: np.ndarray) -> tuple[int, int]:
        return (
            int(round(float(point[0]) * self._scale)),
            int(round(float(point[1]) * self._scale)),
        )

    def add_point(self, point: np.ndarray, *, reuse_existing: bool) -> int:
        key = self._key(point)
        if reuse_existing and key in self._point_index:
            return self._point_index[key]

        index = len(self._points)
        self._points.append(np.asarray(point, dtype=np.float64))
        if reuse_existing:
            self._point_index[key] = index
        return index

    def add_segment(self, start: int, end: int) -> None:
        if start == end:
            return

        canonical = (start, end) if start < end else (end, start)
        if canonical in self._segment_index:
            return

        self._segment_index.add(canonical)
        self._segments.append((start, end))

    def add_loop(self, loop: object, *, reuse_existing: bool) -> np.ndarray:
        clean_loop = _sanitize_loop(loop)
        indices = [self.add_point(point, reuse_existing=reuse_existing) for point in clean_loop]
        for index, next_index in zip(indices, indices[1:] + indices[:1]):
            self.add_segment(index, next_index)
        return clean_loop

    def to_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        points = np.asarray(self._points, dtype=np.float64)
        segments = np.asarray(self._segments, dtype=np.uint32)
        if points.size == 0:
            points = np.empty((0, 2), dtype=np.float64)
        if segments.size == 0:
            segments = np.empty((0, 2), dtype=np.uint32)
        return points, segments


@dataclass(slots=True)
class Domain:
    points: np.ndarray
    segments: np.ndarray | None = None
    holes: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.points = _require_points(self.points, name="points")
        self.segments = _coerce_segments(self.segments)
        self.holes = _coerce_points(self.holes)

    @classmethod
    def from_loops(
        cls,
        outer: object,
        holes: list[object] | tuple[object, ...] | None = None,
        *,
        tolerance: float = 1e-10,
    ) -> Domain:
        builder = _DomainBuilder(tolerance=tolerance)
        builder.add_loop(outer, reuse_existing=False)

        hole_points: list[np.ndarray] = []
        for hole in holes or []:
            hole_loop = builder.add_loop(hole, reuse_existing=False)
            hole_points.append(_loop_interior_point(hole_loop))

        points, segments = builder.to_arrays()
        holes_array = np.asarray(hole_points, dtype=np.float64)
        if holes_array.size == 0:
            holes_array = np.empty((0, 2), dtype=np.float64)
        return cls(points=points, segments=segments, holes=holes_array)


@dataclass(slots=True)
class CoverageGraph:
    points: np.ndarray
    segments: np.ndarray
    region_points: np.ndarray
    region_markers: np.ndarray

    def __post_init__(self) -> None:
        self.points = _require_points(self.points, name="points")
        self.segments = _require_segments(self.segments, name="segments")
        self.region_points = _require_points(self.region_points, name="region_points")
        self.region_markers = _require_markers(self.region_markers, name="region_markers")


@dataclass(slots=True)
class Coverage:
    polygons: tuple[object, ...]
    markers: tuple[int, ...]
    tolerance: float | None = None

    def __init__(
        self,
        polygons: object,
        markers: object,
        *,
        tolerance: float | None = None,
    ) -> None:
        self.polygons = tuple(polygons)
        self.markers = tuple(int(marker) for marker in markers)
        self.tolerance = float(tolerance) if tolerance is not None else None

    def graph(self, *, max_edge_length: float | None = None) -> CoverageGraph:
        points, segments, region_points, region_markers = build_coverage_domain(
            self.polygons,
            self.markers,
            max_edge_length=max_edge_length,
            tolerance=self.tolerance,
        )
        return CoverageGraph(points, segments, region_points, region_markers)


def _domain_mesh(domain: Domain, options: MeshingOptions) -> Mesh:
    raw = _core._generate_raw(
        domain.points,
        domain.segments,
        domain.holes,
        options.min_angle,
        options.max_area,
        options.max_edge_length,
        options.refine,
        options.off_centers,
        options.verbose,
        _acute_mode_value(options.acute_protection),
        options.protect_angle,
        options.max_refine_steps,
        options.max_protection_levels,
    )
    return _mesh_from_raw(raw)


def _coverage_graph_mesh(graph: CoverageGraph, options: MeshingOptions) -> Mesh:
    raw = _core._generate_coverage_raw(
        graph.points,
        graph.segments,
        graph.region_points,
        graph.region_markers,
        options.min_angle,
        options.max_area,
        options.max_edge_length,
        options.refine,
        options.off_centers,
        options.verbose,
        _acute_mode_value(options.acute_protection),
        options.protect_angle,
        options.max_refine_steps,
        options.max_protection_levels,
    )
    return _mesh_from_raw(raw)


def mesh(
    geometry: Domain | CoverageGraph | Coverage,
    *,
    options: MeshingOptions | None = None,
) -> Mesh:
    resolved_options = options or MeshingOptions()

    if isinstance(geometry, Domain):
        return _domain_mesh(geometry, resolved_options)
    if isinstance(geometry, CoverageGraph):
        return _coverage_graph_mesh(geometry, resolved_options)
    if isinstance(geometry, Coverage):
        return _coverage_graph_mesh(
            geometry.graph(max_edge_length=resolved_options.max_edge_length),
            resolved_options,
        )

    raise TypeError(
        "mesh() expects a Domain, CoverageGraph, or Coverage"
    )


def read_domain(path: str | Path | PathLike[str]) -> Domain:
    points, segments, holes = _io.read_domain_file(path)
    return Domain(points=points, segments=segments, holes=holes)
