from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np


def _sanitize_loop(points: object, tol: float = 1e-10) -> np.ndarray:
    if points is None:
        return np.empty((0, 2), dtype=np.float64)

    array = np.asarray(points, dtype=np.float64)
    if array.size == 0:
        return np.empty((0, 2), dtype=np.float64)

    clean: list[np.ndarray] = []
    for point in array:
        if not clean or np.linalg.norm(point - clean[-1]) > tol:
            clean.append(point)

    if len(clean) > 1 and np.linalg.norm(clean[0] - clean[-1]) <= tol:
        clean.pop()

    return np.asarray(clean, dtype=np.float64)


class _CoverageGraphBuilder:
    def __init__(self, tolerance: float = 1e-9):
        self._tolerance = float(tolerance)
        self._scale = 1.0 / tolerance
        self._points: list[np.ndarray] = []
        self._segments: list[tuple[int, int]] = []
        self._point_index: dict[tuple[int, int], list[int]] = {}
        self._segment_index: set[tuple[int, int]] = set()

    def _key(self, point: np.ndarray) -> tuple[int, int]:
        return (
            int(round(float(point[0]) * self._scale)),
            int(round(float(point[1]) * self._scale)),
        )

    def add_point(self, point: np.ndarray, *, reuse_existing: bool = True) -> int:
        key = self._key(point)
        if reuse_existing:
            best_index = None
            best_distance = None
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    neighbor_key = (key[0] + dx, key[1] + dy)
                    for index in self._point_index.get(neighbor_key, []):
                        distance = float(np.linalg.norm(point - self._points[index]))
                        if distance > self._tolerance:
                            continue
                        if best_distance is None or distance < best_distance:
                            best_index = index
                            best_distance = distance
            if best_index is not None:
                return best_index

        index = len(self._points)
        self._points.append(np.asarray(point, dtype=np.float64))
        if reuse_existing:
            self._point_index.setdefault(key, []).append(index)
        return index

    def add_segment(self, start: int, end: int) -> None:
        if start == end:
            return

        canonical = (start, end) if start < end else (end, start)
        if canonical in self._segment_index:
            return

        self._segment_index.add(canonical)
        self._segments.append((start, end))

    def add_loop(self, loop: object) -> None:
        clean_loop = _sanitize_loop(loop)
        if len(clean_loop) < 3:
            return

        indices = [self.add_point(point, reuse_existing=True) for point in clean_loop]
        for index, next_index in zip(indices, indices[1:] + indices[:1]):
            self.add_segment(index, next_index)

    def to_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        points = np.asarray(self._points, dtype=np.float64)
        segments = np.asarray(self._segments, dtype=np.uint32)

        if points.size == 0:
            points = np.empty((0, 2), dtype=np.float64)
        if segments.size == 0:
            segments = np.empty((0, 2), dtype=np.uint32)

        return points, segments


def _iter_lines(geometry: object) -> list[np.ndarray]:
    from shapely.geometry import GeometryCollection, LineString, LinearRing, MultiLineString

    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, (LineString, LinearRing)):
        coords = np.asarray(geometry.coords, dtype=np.float64)
        return [coords] if len(coords) >= 2 else []
    if isinstance(geometry, (GeometryCollection, MultiLineString)):
        lines: list[np.ndarray] = []
        for item in geometry.geoms:
            lines.extend(_iter_lines(item))
        return lines
    return []


def _polygon_interior_seeds(polygon: object, spacing: float | None) -> list[np.ndarray]:
    if spacing is None or spacing <= 0.0:
        return []

    minx, miny, maxx, maxy = polygon.bounds
    if maxx - minx <= spacing or maxy - miny <= spacing:
        return []

    from shapely.geometry import Point
    from shapely.prepared import prep

    prepared = prep(polygon)
    boundary = polygon.boundary
    vertical_step = spacing * math.sqrt(3.0) / 2.0
    boundary_clearance = 0.5 * vertical_step
    seeds: list[np.ndarray] = []

    row = 0
    y = miny + 0.5 * vertical_step
    while y < maxy:
        x_offset = 0.5 * spacing if row % 2 else 0.0
        x = minx + 0.5 * spacing + x_offset
        while x < maxx:
            point = Point(x, y)
            if prepared.contains(point) and boundary.distance(point) >= boundary_clearance * (1.0 - 1e-9):
                seeds.append(np.array([x, y], dtype=np.float64))
            x += spacing
        y += vertical_step
        row += 1

    return seeds


def _resolve_coverage_tolerance(
    polygons: Sequence[object],
    tolerance: float | None,
) -> float:
    if tolerance is not None:
        return float(tolerance)

    bounds = [
        polygon.bounds
        for polygon in polygons
        if polygon is not None and not polygon.is_empty
    ]
    if not bounds:
        return 1e-9

    minx = min(bound[0] for bound in bounds)
    miny = min(bound[1] for bound in bounds)
    maxx = max(bound[2] for bound in bounds)
    maxy = max(bound[3] for bound in bounds)
    diagonal = math.hypot(maxx - minx, maxy - miny)
    return min(max(diagonal * 1e-6, 1e-9), 1e-3)


def build_coverage_domain(
    polygons: Sequence[object],
    markers: Sequence[int],
    *,
    max_edge_length: float | None = None,
    tolerance: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if len(polygons) != len(markers):
        raise ValueError("polygons and markers must have the same length")

    resolved_tolerance = _resolve_coverage_tolerance(polygons, tolerance)
    builder = _CoverageGraphBuilder(tolerance=resolved_tolerance)
    region_points: list[np.ndarray] = []
    region_markers: list[int] = []
    boundary_geometries: list[object] = []

    for polygon, marker in zip(polygons, markers):
        if polygon is None or polygon.is_empty:
            continue

        boundary_geometries.append(polygon.boundary)
        for seed in _polygon_interior_seeds(polygon, max_edge_length):
            builder.add_point(seed, reuse_existing=True)

        region_points.append(
            np.asarray(polygon.representative_point().coords[0], dtype=np.float64)
        )
        region_markers.append(int(marker))

    if boundary_geometries:
        from shapely.ops import unary_union

        noded_boundaries = unary_union(boundary_geometries)
        for line in _iter_lines(noded_boundaries):
            is_closed = (
                len(line) >= 2
                and np.linalg.norm(line[0] - line[-1]) <= resolved_tolerance
            )
            clean_line = _sanitize_loop(line, tol=0.0)
            if len(clean_line) < 2:
                continue
            indices = [builder.add_point(point, reuse_existing=True) for point in clean_line]
            for index, next_index in zip(indices, indices[1:]):
                builder.add_segment(index, next_index)
            if is_closed and len(indices) > 2:
                builder.add_segment(indices[-1], indices[0])

    points, segments = builder.to_arrays()
    region_points_array = np.asarray(region_points, dtype=np.float64)
    region_markers_array = np.asarray(region_markers, dtype=np.int32)

    if region_points_array.size == 0:
        region_points_array = np.empty((0, 2), dtype=np.float64)
    if region_markers_array.size == 0:
        region_markers_array = np.empty((0,), dtype=np.int32)

    return points, segments, region_points_array, region_markers_array
