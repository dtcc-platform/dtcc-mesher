from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from . import _core
from . import _io
from ._plot import plot_mesh

__all__ = [
    "Mesh",
    "QualitySummary",
    "generate",
    "generate_file",
    "plot_mesh",
    "read_domain_file",
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


def generate(
    points: object,
    *,
    segments: object | None = None,
    holes: object | None = None,
    min_angle: float = 20.0,
    max_area: float | None = None,
    refine: bool = True,
    off_centers: bool = False,
    verbose: bool = False,
    acute_protection: Literal["none", "simple", "shell"] = "shell",
    protect_angle: float | None = None,
    max_refine_steps: int = 0,
    max_protection_levels: int = 6,
) -> Mesh:
    if max_area is not None:
        raise NotImplementedError("max_area is not implemented in dtcc_mesher yet")

    points_array = _coerce_points(points)
    segments_array = _coerce_segments(segments)
    holes_array = _coerce_points(holes)

    raw = _core._generate_raw(
        points_array,
        segments_array,
        holes_array,
        min_angle,
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
    )


def read_domain_file(path: str | Path) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    return _io.read_domain_file(path)


def generate_file(path: str | Path, **kwargs) -> Mesh:
    points, segments, holes = read_domain_file(path)
    return generate(points, segments=segments, holes=holes, **kwargs)
