from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from . import Mesh


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def read_domain_file(path: str | Path) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    path = Path(path)
    text = path.read_text(encoding="utf-8").splitlines()
    data_lines = []
    for line in text:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        data_lines.append(stripped)

    if path.suffix == ".pts":
        points = np.array([[float(x), float(y)] for x, y in (line.split() for line in data_lines)], dtype=np.float64)
        return points, None, None

    if path.suffix != ".pslg":
        raise ValueError(f"unsupported input format for {path}")

    cursor = 0
    header = data_lines[cursor].split()
    if len(header) != 2 or header[0] != "vertices":
        raise ValueError(f"invalid PSLG header in {path}")
    num_points = int(header[1])
    cursor += 1

    points = np.array(
        [[float(x), float(y)] for x, y in (data_lines[cursor + i].split() for i in range(num_points))],
        dtype=np.float64,
    )
    cursor += num_points

    header = data_lines[cursor].split()
    if len(header) != 2 or header[0] != "segments":
        raise ValueError(f"missing segments header in {path}")
    num_segments = int(header[1])
    cursor += 1

    segments = np.array(
        [[int(a), int(b)] for a, b in (data_lines[cursor + i].split() for i in range(num_segments))],
        dtype=np.uint32,
    )
    cursor += num_segments

    holes = None
    if cursor < len(data_lines):
        header = data_lines[cursor].split()
        if len(header) != 2 or header[0] != "holes":
            raise ValueError(f"invalid holes header in {path}")
        num_holes = int(header[1])
        cursor += 1
        holes = np.array(
            [[float(x), float(y)] for x, y in (data_lines[cursor + i].split() for i in range(num_holes))],
            dtype=np.float64,
        )

    return points, segments, holes


def write_triangles(mesh: "Mesh", path: str | Path) -> None:
    path = Path(path)
    _ensure_parent(path)
    with path.open("w", encoding="utf-8") as stream:
        for tri in mesh.triangles:
            stream.write(f"{int(tri[0])} {int(tri[1])} {int(tri[2])}\n")


def _triangle_metrics(points: np.ndarray, tri: np.ndarray) -> tuple[float, float, float, float, float]:
    a = points[int(tri[0])]
    b = points[int(tri[1])]
    c = points[int(tri[2])]
    ab = np.linalg.norm(a - b)
    bc = np.linalg.norm(b - c)
    ca = np.linalg.norm(c - a)
    shortest = min(ab, bc, ca)
    longest = max(ab, bc, ca)
    area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) * 0.5
    if area <= 0.0 or shortest <= 0.0:
        raise ValueError("degenerate triangle encountered")

    def angle(left: float, right: float, opposite: float) -> float:
        cosine = (left * left + right * right - opposite * opposite) / (2.0 * left * right)
        cosine = min(1.0, max(-1.0, cosine))
        return float(np.degrees(np.arccos(cosine)))

    angle_a = angle(ab, ca, bc)
    angle_b = angle(ab, bc, ca)
    angle_c = angle(bc, ca, ab)
    min_angle = min(angle_a, angle_b, angle_c)
    max_angle = max(angle_a, angle_b, angle_c)
    edge_ratio = longest / shortest
    circumradius = (ab * bc * ca) / (4.0 * area)
    radius_edge_ratio = circumradius / shortest
    return area, min_angle, max_angle, edge_ratio, radius_edge_ratio


def write_quality_csv(mesh: "Mesh", path: str | Path) -> None:
    path = Path(path)
    _ensure_parent(path)
    with path.open("w", encoding="utf-8") as stream:
        stream.write("tri_id,v0,v1,v2,area,min_angle_deg,max_angle_deg,edge_ratio,radius_edge_ratio\n")
        for i, tri in enumerate(mesh.triangles):
            area, min_angle, max_angle, edge_ratio, radius_edge_ratio = _triangle_metrics(mesh.points, tri)
            stream.write(
                f"{i},{int(tri[0])},{int(tri[1])},{int(tri[2])},{area:.17g},{min_angle:.17g},{max_angle:.17g},{edge_ratio:.17g},{radius_edge_ratio:.17g}\n"
            )


def write_quality_summary(mesh: "Mesh", path: str | Path) -> None:
    path = Path(path)
    _ensure_parent(path)
    summary = mesh.summary
    with path.open("w", encoding="utf-8") as stream:
        for key in (
            "triangle_count",
            "point_count",
            "input_point_count",
            "steiner_point_count",
            "segment_split_point_count",
            "triangle_split_point_count",
            "protected_corner_count",
            "exempt_triangle_count",
            "area_min",
            "area_mean",
            "area_max",
            "min_angle_deg_min",
            "min_angle_deg_mean",
            "min_angle_deg_max",
            "edge_ratio_min",
            "edge_ratio_mean",
            "edge_ratio_max",
            "radius_edge_ratio_min",
            "radius_edge_ratio_mean",
            "radius_edge_ratio_max",
            "count_min_angle_lt_20",
            "count_min_angle_lt_30",
        ):
            stream.write(f"{key}={getattr(summary, key)}\n")


def write_svg(mesh: "Mesh", path: str | Path) -> None:
    path = Path(path)
    _ensure_parent(path)

    min_xy = mesh.points.min(axis=0)
    max_xy = mesh.points.max(axis=0)
    span = np.maximum(max_xy - min_xy, 1.0)
    margin = 24.0
    scale = 560.0 / float(max(span[0], span[1]))
    canvas_w = float(span[0] * scale + 2.0 * margin)
    canvas_h = float(span[1] * scale + 2.0 * margin)

    def map_point(point: np.ndarray) -> tuple[float, float]:
        x = margin + float(point[0] - min_xy[0]) * scale
        y = canvas_h - margin - float(point[1] - min_xy[1]) * scale
        return x, y

    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{canvas_w:.0f}\" height=\"{canvas_h:.0f}\" viewBox=\"0 0 {canvas_w:.0f} {canvas_h:.0f}\">\n"
        )
        stream.write("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n")

        for tri in mesh.triangles:
            a = map_point(mesh.points[int(tri[0])])
            b = map_point(mesh.points[int(tri[1])])
            c = map_point(mesh.points[int(tri[2])])
            stream.write(
                f"  <polygon points=\"{a[0]:.6f},{a[1]:.6f} {b[0]:.6f},{b[1]:.6f} {c[0]:.6f},{c[1]:.6f}\" fill=\"none\" stroke=\"black\" stroke-width=\"1.25\"/>\n"
            )

        for seg in mesh.segments:
            a = map_point(mesh.points[int(seg[0])])
            b = map_point(mesh.points[int(seg[1])])
            stream.write(
                f"  <line x1=\"{a[0]:.6f}\" y1=\"{a[1]:.6f}\" x2=\"{b[0]:.6f}\" y2=\"{b[1]:.6f}\" stroke=\"#b22222\" stroke-width=\"3\"/>\n"
            )

        if len(mesh.points) <= 16:
            for i, point in enumerate(mesh.points):
                px, py = map_point(point)
                stream.write(
                    f"  <text x=\"{px + 5.0:.6f}\" y=\"{py - 5.0:.6f}\" font-family=\"monospace\" font-size=\"12\" fill=\"black\">{i}</text>\n"
                )

        stream.write("</svg>\n")
