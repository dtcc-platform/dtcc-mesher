from __future__ import annotations

import dtcc_mesher as dm

points = [
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
    (0.5, 0.5),
]

mesh = dm.generate(points)
print(f"points={mesh.summary.point_count} triangles={mesh.summary.triangle_count}")
print(f"min_angle_deg_min={mesh.summary.min_angle_deg_min:.2f}")
