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

## CLI

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

mesh = dm.generate(
    [
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (0.5, 0.5),
    ]
)

print(mesh.points.shape)
print(mesh.triangles.shape)
print(mesh.summary.min_angle_deg_min)
mesh.write_svg("mesh.svg")
```

```python
import dtcc_mesher as dm

mesh = dm.generate_file(
    "tests/cases/square_hole_domain.pslg",
    min_angle=25.0,
    max_area=0.5,
    off_centers=True,
)
mesh.write_quality_summary("square_hole.summary.txt")
mesh.show(title="square_hole_domain")
```

```sh
python -m dtcc_mesher tests/cases/square_hole_domain.pslg build/out/square_hole_py
```

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
