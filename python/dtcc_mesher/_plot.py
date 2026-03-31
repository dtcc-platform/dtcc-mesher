from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import Mesh


def plot_mesh(mesh: "Mesh", ax=None, show_segments: bool = True):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for Mesh.plot(); install dtcc_mesher[plot]") from exc

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

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    if created:
        return ax.figure, ax
    return ax
