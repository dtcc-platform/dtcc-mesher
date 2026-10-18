"""Constraint recovery must try legal crossed edges before adding points."""
import json
from pathlib import Path
import numpy as np
import dtcc_mesher as dm


def test_concave_crossing_quadrilateral_does_not_force_recovery_points():
    data=json.loads((Path(__file__).parents[1]/'cases/recovery_concave_crossing.json').read_text())
    domain=dm.Domain.from_loops(np.array(data['outer']),[np.array(p) for p in data['holes']])
    result=dm.mesh(domain,options=dm.MeshingOptions(refine=False,acute_protection='none'))
    assert np.array_equal(result.points,domain.points)
    assert len(result.triangles)==len(domain.points)+2*len(data['holes'])-2
    assert {tuple(sorted(e)) for e in result.segments}=={tuple(sorted(e)) for e in domain.segments}
    assert set(np.unique(result.triangles))==set(range(len(domain.points)))
    assert dm.improve(result,max_added_points=0).mesh.points.shape==domain.points.shape
