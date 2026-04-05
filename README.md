# DTCC Mesher

MIT Licensed Two-Dimensional Quality Mesh Generator.

## Features

- point-set Delaunay triangulation from `.pts` files
- constrained meshing for `.pslg` domains
- quality refinement with midpoint segment splitting and Steiner insertion
- optional off-centers
- acute-corner protection
- SVG, CSV, and text reports
- C library
- native CLI
- Python package with NumPy arrays
- CMake package export and `pkg-config` metadata

## Limitations

- `.pts` inputs are not refined
- refinement applies to `.pslg` inputs

## Build

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix ./install
```

Developer shortcuts:

```sh
make
make test
```

## Python Install

```sh
python -m pip install .
python -m pip install ".[plot]"
python -m pip install ".[dev]"
```

## Native CLI

The native CLI is the C executable built by CMake:

```sh
./build/dtcc_mesher tests/cases/square_center5.pts build/out/square_center5
./build/dtcc_mesher tests/cases/square_hole_domain.pslg build/out/square_hole
./build/dtcc_mesher --off-centers tests/cases/city_tight_downtown_domain.pslg build/out/city_tight
```

Interface:

```text
dtcc_mesher [options] input.(pts|pslg) outbase
```

Options:

- `-h`, `--help`
- `--version`
- `-v`, `--verbose`
- `--refine`, `--no-refine`
- `--off-centers`, `--no-off-centers`
- `--acute-protection`
- `--simple-acute-protection`
- `--shell-acute-protection`
- `--no-acute-protection`
- `--min-angle deg`
- `--max-area area`
- `--protect-angle deg`
- `--max-refine-steps n`
- `--max-protection-levels n`

Outputs:

- `outbase.tri`
- `outbase.svg`
- `outbase.metrics.csv`
- `outbase.summary.txt`

## C API

Public headers:

- `dtcc_mesher/dtcc_mesher.h`
- `dtcc_mesher/dtcc_mesher_io.h`
- `dtcc_mesher/dtcc_mesher_version.h`

```c
#include <stdio.h>
#include "dtcc_mesher/dtcc_mesher.h"

int main(void) {
    dtcc_mesher_point points[] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    dtcc_mesher_domain domain = {
        .points = points,
        .num_points = 5,
        .segments = NULL,
        .num_segments = 0,
        .holes = NULL,
        .num_holes = 0
    };
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_error error;

    dtcc_mesher_options_init(&options);
    if (dtcc_mesher_generate(&domain, &options, &mesh, &error) != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "%s\n", error.message);
        return 1;
    }

    dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    printf("triangles=%zu min_angle=%.2f\n", summary.triangle_count, summary.min_angle_deg_min);
    dtcc_mesher_mesh_free(&mesh);
    return 0;
}
```

## Python

```python
import dtcc_mesher as dm

mesh = dm.mesh(
    dm.Domain(
        points=[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, 0.5),
        ]
    )
)

print(mesh.points.shape)
print(mesh.triangles.shape)
print(mesh.summary.min_angle_deg_min)
mesh.write_svg("mesh.svg")
```

```python
import dtcc_mesher as dm

domain = dm.read_domain("tests/cases/square_hole_domain.pslg")
mesh = dm.mesh(
    domain,
    options=dm.MeshingOptions(
        min_angle=25.0,
        max_edge_length=0.75,
        off_centers=True,
    ),
)
mesh.write_quality_summary("square_hole.summary.txt")
mesh.show(title="square_hole_domain")
```

```python
import dtcc_mesher as dm
from shapely.geometry import box

mesh = dm.mesh(
    dm.Coverage(
        [box(0.0, 0.0, 0.5, 1.0), box(0.5, 0.0, 1.0, 1.0)],
        markers=[10, 20],
    ),
    options=dm.MeshingOptions(max_edge_length=0.25, refine=False),
)
```

Public Python API:

- `dm.mesh(...)` is the single entry point for meshing
- `dm.Domain(...)` and `dm.Domain.from_loops(...)` are for PSLG / point-segment-hole inputs
- `dm.Coverage(...)` is for polygon coverage inputs with region markers
- `dm.CoverageGraph(...)` is the low-level coverage graph input when you already have noded segments and region markers
- `dm.read_domain(...)` reads `.pts` / `.pslg` files into a `Domain`
- `dm.MeshingOptions(...)` carries the meshing controls
- `dm.Mesh` carries `points`, `triangles`, `segments`, `summary`, and optional `markers`
- `dm.plot_mesh(...)`, `dm.plot_mesh_with_summary(...)`, and `dm.show_mesh(...)` are plotting helpers

Recommended usage pattern:

1. Build or read a `Domain` / `Coverage`
2. Create `MeshingOptions` only when you want non-default controls
3. Call `mesh(...)`

`MeshingOptions` fields:

- `min_angle`
- `max_edge_length`:
  preferred global size control, interpreted as the maximum target constrained-edge
  length / nominal mesh spacing in the input plane. Use `None` for unrestricted size.
- `max_area`:
  lower-level area cap for backends that refine by area; most callers should prefer
  `max_edge_length`
- `refine`
- `off_centers`
- `verbose`
- `acute_protection`
- `protect_angle`
- `max_refine_steps`
- `max_protection_levels`

## Python CLI

The Python package also exposes a lightweight CLI:

```sh
python -m dtcc_mesher tests/cases/square_hole_domain.pslg build/out/square_hole_py
```

Options:

- `-h`, `--help`
- `--min-angle deg`
- `--max-area area`
- `--max-edge-length length`
- `--no-refine`
- `--off-centers`
- `--acute-protection {none,simple,shell}`

Plotting requires:

```sh
python -m pip install ".[plot]"
```

## Demos

- `demos/c/basic_generate.c`
- `demos/c/quality_report.c`
- `demos/python/basic_generate.py`
- `demos/python/plot_mesh.py`
- `demos/cli/README.md`

Interactive plotting demo:

```sh
python demos/python/plot_mesh.py
python demos/python/plot_mesh.py --list
python demos/python/plot_mesh.py city_tight_downtown_domain --off-centers
```

## Tests

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

```sh
python -m pip install ".[dev]"
python -m pytest tests/python
```

## Layout

```text
include/dtcc_mesher/   public C headers
src/api/               public C API wrappers
src/core/              internal mesh generator sources
src/third_party/       public-domain predicates
cli/                   native CLI
python/                pybind11 binding and Python package
demos/                 C, Python, and CLI demos
tests/                 algorithmic, API, CLI, and Python tests
```

## Third-Party Code

DTCC Mesher is MIT licensed. `src/third_party/predicates.c` is copied verbatim from Jonathan Richard Shewchuk's public-domain robust predicates implementation. All other project files are part of DTCC Mesher.
