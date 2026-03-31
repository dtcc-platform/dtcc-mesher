from __future__ import annotations

import dtcc_mesher as dm

mesh = dm.generate_file("tests/cases/square_hole_domain.pslg")
mesh.write_svg("demo_plot_mesh.svg")
print("wrote demo_plot_mesh.svg")
