from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import Mesh, QualitySummary


def _load_pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting; install dtcc_mesher[plot]") from exc

    return plt


def format_summary_lines(summary: "QualitySummary") -> list[str]:
    return [
        f"triangles: {summary.triangle_count}",
        f"points: {summary.point_count}",
        f"input points: {summary.input_point_count}",
        f"steiner points: {summary.steiner_point_count}",
        f"segment splits: {summary.segment_split_point_count}",
        f"triangle splits: {summary.triangle_split_point_count}",
        f"protected corners: {summary.protected_corner_count}",
        f"exempt triangles: {summary.exempt_triangle_count}",
        f"area min/mean/max: {summary.area_min:.6g} / {summary.area_mean:.6g} / {summary.area_max:.6g}",
        (
            "min angle min/mean/max: "
            f"{summary.min_angle_deg_min:.4f} / {summary.min_angle_deg_mean:.4f} / {summary.min_angle_deg_max:.4f}"
        ),
        f"edge ratio min/mean/max: {summary.edge_ratio_min:.6g} / {summary.edge_ratio_mean:.6g} / {summary.edge_ratio_max:.6g}",
        (
            "radius-edge ratio min/mean/max: "
            f"{summary.radius_edge_ratio_min:.6g} / {summary.radius_edge_ratio_mean:.6g} / {summary.radius_edge_ratio_max:.6g}"
        ),
        f"count min angle < 20: {summary.count_min_angle_lt_20}",
        f"count min angle < 30: {summary.count_min_angle_lt_30}",
    ]


def plot_mesh(
    mesh: "Mesh",
    ax=None,
    show_segments: bool = True,
    show_points: bool = False,
    title: str | None = None,
):
    plt = _load_pyplot()
    created = False
    if ax is None:
        _, ax = plt.subplots()
        created = True

    for tri in mesh.triangles:
        coords = mesh.points[tri]
        cycle = coords[[0, 1, 2, 0]]
        ax.plot(cycle[:, 0], cycle[:, 1], color="black", linewidth=0.8)

    if show_segments and len(mesh.segments) != 0:
        for seg in mesh.segments:
            coords = mesh.points[seg]
            ax.plot(coords[:, 0], coords[:, 1], color="#b22222", linewidth=1.8)

    if show_points:
        ax.scatter(mesh.points[:, 0], mesh.points[:, 1], color="#1f4e79", s=18.0, zorder=3)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.margins(0.05)
    if title is not None:
        ax.set_title(title)

    if created:
        return ax.figure, ax
    return ax


def plot_mesh_with_summary(
    mesh: "Mesh",
    *,
    title: str | None = None,
    show_segments: bool = True,
    show_points: bool = False,
    figsize: tuple[float, float] = (12.0, 7.0),
):
    plt = _load_pyplot()
    fig, (ax_mesh, ax_summary) = plt.subplots(
        1,
        2,
        figsize=figsize,
        gridspec_kw={"width_ratios": [3.0, 1.6]},
        constrained_layout=True,
    )

    plot_mesh(mesh, ax=ax_mesh, show_segments=show_segments, show_points=show_points, title=title)

    ax_summary.axis("off")
    ax_summary.set_title("Summary", loc="left")
    ax_summary.text(
        0.0,
        1.0,
        "\n".join(format_summary_lines(mesh.summary)),
        va="top",
        ha="left",
        family="monospace",
        fontsize=10,
    )

    return fig, (ax_mesh, ax_summary)


def show_mesh(
    mesh: "Mesh",
    *,
    title: str | None = None,
    show_segments: bool = True,
    show_points: bool = False,
    block: bool = True,
):
    plt = _load_pyplot()
    fig, axes = plot_mesh_with_summary(
        mesh,
        title=title,
        show_segments=show_segments,
        show_points=show_points,
    )
    plt.show(block=block)
    return fig, axes
