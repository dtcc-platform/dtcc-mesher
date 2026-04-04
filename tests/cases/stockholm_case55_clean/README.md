Stockholm case `55` exported from `dtcc-core` after footprint cleaning.

What is stored here:
- One `.pslg` file per actual `dtcc_mesher.mesh(dtcc_mesher.Domain(...))` call used by the flat-mesh path for this tile.
- `groundXX_domain.pslg` files for the ground-domain components after subtracting cleaned buildings from the tile bounds.
- `buildingXX_domain.pslg` files for the cleaned building domains.
- `manifest.json` with tile metadata, mesher settings, per-domain sizes, and the current `dtcc_mesher` success/failure result for each domain.

Source configuration:
- Origin: `dtcc-core/sandbox/mesh_quality_survey.py`
- Tile: case `55`
- Bounds: `(675000, 6581000) -> (675500, 6581500)` in `EPSG:3006`
- Cleaning path: `build_city_flat_mesh(..., lod=LOD0, min_building_detail=0.5, min_building_area=15.0, merge_tolerance=0.5, merge_buildings=True)`
- Mesher input settings: `max_mesh_size=10.0`, `min_angle=25.0`, `refine=False`

This corpus is intentionally a current stress case. At export time, `dtcc_mesher` failed on 4 of the 89 stored domains.
