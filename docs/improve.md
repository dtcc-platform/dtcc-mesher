# Stateless bounded shape improvement

`mesh()` remains constrained Delaunay generation/refinement. The separate
`improve()` operation optimizes an existing valid constrained planar mesh and
**does not promise Delaunay output**. It never moves or removes input vertices.

```python
import dtcc_mesher as dm
mesh = dm.mesh(domain, options=dm.MeshingOptions(refine=False,
                                               acute_protection="none"))
zero = dm.improve(mesh, max_added_points=0)
trial = dm.improve(zero.mesh, max_added_points=1,
                   splittable_segments=[segment_index], allow_interior=False)
if trial.insertion_kind == "segment":
    parent = trial.parent_segment
    t = trial.segment_parameter  # original segment a + t*(b-a)
```

`max_added_points` accepts only integers 0 or 1. Boundary and internal constrained
segments remain immutable unless their input index is explicitly permitted.
`allow_interior=True` permits an interior trial; `False` restricts insertion to
permitted segments. This one flag lets a caller enforce a separate interior cap
and evaluate shared-edge costs without an interior winner hiding an edge trial.
No permission means no segment insertion. Zero points always permits legal flips.

The result contains an ordinary owned `Mesh`, `insertion_kind` (`none`, `interior`,
`segment`), an appended `inserted_point` index or None, parent segment index and
oriented parameter, and one `segment_parents` entry per output segment. Child
segments preserve their parent's direction. All original point indices and
coordinates are retained exactly. No implicit recovery or refinement points are
allowed. Counts in the result mesh treat all supplied vertices as this operation's
input; protection/exemption metadata is not a constraint. Triangle semantic
markers are unsupported by this operation: split semantic domains before calling.

Diagnostics contain before/after shape scores and existing quality summaries,
accepted flip count, candidate evaluations and termination. `no_beneficial_trial`
is a normal result, not a geometry-impossibility proof. Invalid input raises
ValueError; an unexpected native numerical/topology failure raises RuntimeError;
allocation failure raises MemoryError. Input arrays are never mutated.

`segment_trials` contains at most four beneficial, validated segment trials,
ranked by native gain per added triangle (including the segment winner). Each
entry has `parent_segment` and `segment_parameter`; no additional mesh is stored.
An application coordinating incident domains can retessellate them at one of
these parameters and choose by combined cost/quality. A local winner can harm a
neighbour while another trial passes. Zero-point calls return an empty tuple.
This remains a zero-or-one-point operation; the records do not commit points.

## Bounds and criterion

The finite search considers at most 128 distinct candidates, informed by the 12
worst-q triangles. Interior candidates are centroids and equilateral off-centres.
Each allowed segment first receives a midpoint and projections from at most six
worst-q triangles in its endpoint star, with stable triangle-index ties. Global
poor-triangle projections follow. This preserves useful local proposals when
remote acute triangles consume the global selection. Centroids and directly
opposite vertices inform positions; there is no uniform boundary grid. Proposals use the
native numerical coincidence policy, with translated coordinates internally.
Each candidate is applied to a private native mesh using the existing insertion,
segment-registry, legal-flip and validation routines.

`q = 4 sqrt(3) area / sum(edge²)`, and shape-defect energy is
`E = sum((1-q)^2)`. The public score's `penalty` field reports E. This replaces
the unreleased severe-only `sum(max(0, .25-q))` meaning, which cannot reward
moderate fans when all q exceed .25. E is zero for equilateral triangles and
counts the defect of every resulting element, so additional triangles have a
cost even before the explicit gain threshold. Insertion must reduce E by more than .015 per added
triangle. It must not lower the minimum q of the changed region or increase the
absolute count or area below 5° or 10°. An unchanged remote worst triangle is
allowed. This guards moving defects while allowing local progress. Legal flips
increase the worse q of their two triangles without increasing E or these severe
defects. A pass accepts at most `min(8*T+1, 4096)` flips. No minimum-angle-to-
completion loop, size cap, protection shells or cascaded insertion runs.

Scope is a conforming positive-oriented planar embedding with all naked edges
listed as segments, at most 4096 points and 4096 triangles. Holes, fixed internal
constraints and disconnected components are supported. Import checks indices,
finite coordinates, duplicate vertices/triangles, orientation, edge pairing,
constraint presence, crossings, overlap and T junctions. Unsupported input fails
explicitly; input is not silently repaired. After native local changes, topology,
constraints and area are checked without requiring Delaunay edges.

## C interface and ownership

Include `dtcc_mesher/dtcc_mesher_improve.h`. The additive entry point is:

```c
dtcc_mesher_status dtcc_mesher_improve(
    const dtcc_mesher_mesh *input, unsigned int max_added_points,
    int allow_interior, const size_t *splittable_segments, size_t count,
    dtcc_mesher_improvement *result, dtcc_mesher_error *error);
```

Initialize an empty result with `dtcc_mesher_improvement_init`; free all returned
allocations, including unchanged copies, with `dtcc_mesher_improvement_free`.
Input is borrowed for the call and must not alias output storage. On failure the
result is empty (invalid aliasing arguments are rejected without touching storage).
`SIZE_MAX` represents absent point/parent indices in C. The Python binding copies
native arrays into NumPy-owned arrays, frees the native result on every exit, and
releases the GIL during native computation.

The public additions are this operation, its init/free helpers, insertion-kind
enum, shape-score, segment-trial and improvement-result structs, and Python `improve` and
`ImprovementResult`. The bounded segment-trial array/count extends the unreleased
improvement result; C consumers must rebuild with the matching header and library
because its struct size changes. It adds no allocations or operation arguments.
Existing generation options/structs, ordinary generation APIs and
CDT validators have no changes for this capability. Internal export and segment-
registry helpers are reused; mutable topology is not exposed publicly.
