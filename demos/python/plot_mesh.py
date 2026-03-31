from __future__ import annotations

import argparse
from pathlib import Path

import dtcc_mesher as dm

REPO_ROOT = Path(__file__).resolve().parents[2]
CASES_DIR = REPO_ROOT / "tests" / "cases"
INVALID_CASES = {
    "collinear4",
    "duplicate5",
    "intersecting_segments",
}


def available_cases() -> list[Path]:
    paths = sorted(CASES_DIR.glob("*.pts")) + sorted(CASES_DIR.glob("*.pslg"))
    return [path for path in sorted(paths, key=lambda item: item.name) if path_stem(path) not in INVALID_CASES]


def path_stem(path: Path) -> str:
    return path.name.rsplit(".", 1)[0]


def resolve_case(selection: str | None, cases: list[Path]) -> Path:
    by_name = {path.name: path for path in cases}
    by_stem = {path_stem(path): path for path in cases}

    if selection is None:
        print("Available cases:")
        for index, path in enumerate(cases, start=1):
            print(f"  {index:2d}. {path.name}")
        raw = input("Select a case by number or name: ").strip()
        if raw.isdigit():
            index = int(raw)
            if 1 <= index <= len(cases):
                return cases[index - 1]
        selection = raw

    if selection in by_name:
        return by_name[selection]
    if selection in by_stem:
        return by_stem[selection]

    raise SystemExit(f"unknown case: {selection}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a bundled DTCC Mesher case and show it with matplotlib")
    parser.add_argument("case", nargs="?", help="Case name, file name, or menu selection")
    parser.add_argument("--list", action="store_true", help="List available runnable cases and exit")
    parser.add_argument("--show-points", action="store_true", help="Draw point markers")
    parser.add_argument("--off-centers", action="store_true", help="Use off-centers")
    parser.add_argument("--no-refine", action="store_true", help="Disable refinement")
    parser.add_argument(
        "--acute-protection",
        choices=("none", "simple", "shell"),
        default="shell",
        help="Acute-corner protection mode",
    )
    parser.add_argument("--min-angle", type=float, default=20.0, help="Minimum angle target")
    parser.add_argument("--quiet", action="store_true", help="Disable mesher log output")
    return parser.parse_args()


def format_summary_lines(summary) -> list[str]:
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


def show_mesh(mesh, *, title: str, show_points: bool) -> None:
    import matplotlib.pyplot as plt

    fig, (ax_mesh, ax_summary) = plt.subplots(
        1,
        2,
        figsize=(12.0, 7.0),
        gridspec_kw={"width_ratios": [3.0, 1.6]},
        constrained_layout=True,
    )

    mesh.plot(ax=ax_mesh)
    if show_points:
        ax_mesh.scatter(mesh.points[:, 0], mesh.points[:, 1], color="#1f4e79", s=18.0, zorder=3)
    ax_mesh.set_title(title)
    ax_mesh.margins(0.05)

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

    if "agg" in plt.get_backend().lower():
        fig.canvas.draw()
        return

    plt.show()


def main() -> int:
    args = parse_args()
    cases = available_cases()

    if args.list:
        for path in cases:
            print(path.name)
        return 0

    case_path = resolve_case(args.case, cases)

    print(f"case: {case_path.name}")
    print(f"path: {case_path}")
    print()

    mesh = dm.generate_file(
        case_path,
        min_angle=args.min_angle,
        refine=not args.no_refine,
        off_centers=args.off_centers,
        acute_protection=args.acute_protection,
        verbose=not args.quiet,
    )

    print("summary:")
    for line in format_summary_lines(mesh.summary):
        print(f"  {line}")

    title = f"{case_path.name}  |  triangles={mesh.summary.triangle_count}  points={mesh.summary.point_count}"
    show_mesh(mesh, title=title, show_points=args.show_points)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
