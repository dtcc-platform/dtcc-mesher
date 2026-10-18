"""Stable acute protection for coordinated remeshing of fixed domains."""
import numpy as np
import dtcc_mesher as dm


def test_corner_reuse_does_not_recursively_shrink_shell():
    d=dm.Domain.from_loops(np.array([[0.,0.],[.25,0.],[0.,.006]]))
    a=dm.mesh(d,options=dm.MeshingOptions(min_angle=20,max_edge_length=8,max_area=32))
    b=dm.mesh(dm.Domain(a.points[np.unique(a.segments)],np.searchsorted(np.unique(a.segments),a.segments)),
              options=dm.MeshingOptions(min_angle=20,max_edge_length=8,max_area=32,reuse_acute_corners=True))
    # Reusing a shell must not restart its six dyadic protective levels.
    def shortest(m):return np.linalg.norm(m.points[m.segments[:,1]]-m.points[m.segments[:,0]],axis=1).min()
    assert shortest(b)>=shortest(a)*.99
    assert b.quality_exempt.sum()==b.summary.exempt_triangle_count
    assert b.summary.protected_corner_count>0
    np.testing.assert_allclose(b.summary.area_mean*b.summary.triangle_count,.00075,rtol=1e-10)


def test_reuse_unsplit_acute_input_is_bounded_and_explicit():
    d=dm.Domain.from_loops(np.array([[0.,0.],[1.,0.],[0.,.01]]))
    a=dm.mesh(d,options=dm.MeshingOptions(reuse_acute_corners=True,max_refine_steps=1000))
    assert a.summary.protected_corner_count>0
    assert a.quality_exempt.sum()==a.summary.exempt_triangle_count


def test_default_options_keep_corner_reuse_disabled():
    assert not dm.MeshingOptions().reuse_acute_corners
