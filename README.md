# DTCC Mesher

DTCC Mesher provides a free/open-source two-dimensional quality mesh
generator. The purpose is to provide constrained Delaunay
triangulation and quality refinement capabilities similar to Triangle,
while relying on an independent implementation based on published
algorithms and reports rather than on Triangle source code.

This project is part of the
[Digital Twin Platform (DTCC Platform)](https://github.com/dtcc-platform/)
developed at the
[Digital Twin Cities Centre](https://dtcc.chalmers.se/)
supported by Sweden's Innovation Agency Vinnova under Grant No. 2019-421 00041.

![City mesh](docs/images/city_mesh.png)

## Documentation

The public interfaces are:

- the Python package in `python/dtcc_mesher/`
- the public C headers in `include/dtcc_mesher/`
- the native CLI built from `cli/main.c`

The main Python entry points are `mesh(...)`, `Domain`, `Coverage`, `CoverageGraph`,
`MeshingOptions`, and `Mesh`.

## Features

- point-set Delaunay triangulation from `.pts` files
- constrained meshing for `.pslg`, `Coverage`, and `CoverageGraph` inputs
- standard constrained-Delaunay refinement with midpoint segment splitting and
  circumcenter Steiner insertion
- acute-corner protection for small PSLG angles
- C library, native CLI, and Python package

## Limitations

- `.pts` inputs are triangulated but not refined

## Installation

### Python

```sh
python -m pip install .
python -m pip install ".[plot]"
python -m pip install ".[dev]"
```

### Native build

```sh
cmake -S . -B build
cmake --build build
cmake --install build --prefix ./install
```

### Makefile shortcuts

- `make` - configure and build
- `make test` - run the test suite

## Basic usage

Native CLI:

```sh
./build/dtcc_mesher tests/cases/square_hole_domain.pslg build/out/square_hole
```

Python:

```python
import dtcc_mesher as dm

domain = dm.read_domain("tests/cases/square_hole_domain.pslg")
mesh = dm.mesh(domain, options=dm.MeshingOptions(min_angle=25.0, max_edge_length=0.75))
mesh.write_quality_summary("square_hole.summary.txt")
```

For plotting support, install `".[plot]"`. For the complete CLI flags, run
`./build/dtcc_mesher --help` or `python -m dtcc_mesher --help`.

## Algorithmic basis and provenance

DTCC Mesher does not vendor Triangle meshing source code. The only copied third-party
source file in this repository is `src/third_party/predicates.c`, which is Jonathan
Richard Shewchuk's public-domain robust geometric predicates implementation. The mesher
core is instead guided by published work on incremental Delaunay triangulation,
constrained Delaunay refinement, robust predicates, and mesh generation for domains
with small angles.

## References

- Adrian Bowyer, [Computing Dirichlet Tessellations](https://doi.org/10.1093/comjnl/24.2.162), 1981.
- David F. Watson, [Computing the n-dimensional Delaunay Tessellation with Application to Voronoi Polytopes](https://doi.org/10.1093/comjnl/24.2.167), 1981.
- Jim Ruppert, [A Delaunay Refinement Algorithm for Quality 2-Dimensional Mesh Generation](https://doi.org/10.1006/jagm.1995.1021), 1995.
- Jonathan Richard Shewchuk, [Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator](https://www.cs.cmu.edu/~quake/tripaper/triangle0.html), 1996.
- Jonathan Richard Shewchuk, [Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates](https://www.cs.cmu.edu/~quake/robust.html), 1997.
- Jonathan Richard Shewchuk, [Delaunay Refinement Mesh Generation](https://www.cs.cmu.edu/~quake-papers/delaunay-refinement.pdf), 1997.
- Jonathan Richard Shewchuk, [Mesh Generation for Domains with Small Angles](https://people.eecs.berkeley.edu/~jrs/papers/angle.pdf), 2000.
- Jonathan Richard Shewchuk, [Delaunay Refinement Algorithms for Triangular Mesh Generation](https://doi.org/10.1016/S0925-7721(01)00047-5), 2002.

## Authors (in order of appearance)

* [Anders Logg](http://anders.logg.org), with the help of Codex et al.

## License

This project is licensed under the
[MIT license](https://opensource.org/licenses/MIT).

## Community guidelines

Comments, contributions, and questions are welcome. Please engage with
us through Issues, Pull Requests, and Discussions on our GitHub page.
