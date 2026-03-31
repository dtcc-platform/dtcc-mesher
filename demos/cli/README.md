# CLI Demos

These commands assume you already built the project with CMake:

```sh
cmake -S . -B build
cmake --build build
```

Generate a simple point-set triangulation:

```sh
./build/dtcc_mesher tests/cases/square_center5.pts demos/out/square_center5
```

Generate a PSLG mesh with quality refinement:

```sh
mkdir -p demos/out
./build/dtcc_mesher tests/cases/square_hole_domain.pslg demos/out/square_hole
```

Try the tougher downtown stress case:

```sh
./build/dtcc_mesher tests/cases/city_tight_downtown_domain.pslg demos/out/city_tight_downtown
```

Each run writes:

- `*.tri`
- `*.svg`
- `*.metrics.csv`
- `*.summary.txt`
