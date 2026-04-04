from __future__ import annotations

import argparse

from . import generate_file


def main() -> int:
    parser = argparse.ArgumentParser(description="DTCC Mesher Python CLI")
    parser.add_argument("input_path", help="Path to a .pts or .pslg file")
    parser.add_argument("out_base", help="Output file prefix")
    parser.add_argument("--min-angle", type=float, default=20.0, help="Target minimum angle in degrees")
    parser.add_argument("--max-area", type=float, default=None, help="Maximum triangle area for PSLG refinement")
    parser.add_argument("--max-edge-length", type=float, default=None, help="Maximum PSLG segment length before meshing")
    parser.add_argument("--no-refine", action="store_true", help="Disable PSLG refinement")
    parser.add_argument("--off-centers", action="store_true", help="Use off-centers for bad-triangle insertion")
    parser.add_argument(
        "--acute-protection",
        choices=("none", "simple", "shell"),
        default="shell",
        help="Acute-corner protection mode",
    )
    args = parser.parse_args()

    mesh = generate_file(
        args.input_path,
        min_angle=args.min_angle,
        max_area=args.max_area,
        max_edge_length=args.max_edge_length,
        refine=not args.no_refine,
        off_centers=args.off_centers,
        acute_protection=args.acute_protection,
    )
    mesh.write_triangles(f"{args.out_base}.tri")
    mesh.write_svg(f"{args.out_base}.svg")
    mesh.write_quality_csv(f"{args.out_base}.metrics.csv")
    mesh.write_quality_summary(f"{args.out_base}.summary.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
