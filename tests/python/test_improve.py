"""Explicit bounded shape optimization independent of ordinary CDT generation."""
import numpy as np
import pytest
import dtcc_mesher as dm


def mesh(outer, holes=None):
    return dm.mesh(dm.Domain.from_loops(np.asarray(outer,float),holes),
                   options=dm.MeshingOptions(refine=False,acute_protection='none'))


def fan(width=20.,height=8.):
    m=mesh([[x,0] for x in np.linspace(0,width,21)]+[[width,height],[0,height]])
    top=int(np.flatnonzero(np.all(m.points[m.segments,1]==height,axis=1))[0])
    return m,top


def area(m):
    p=m.points[m.triangles];u=p[:,1]-p[:,0];v=p[:,2]-p[:,0]
    return abs(u[:,0]*v[:,1]-u[:,1]*v[:,0]).sum()/2


def test_rectangle_zero_points_input_immutable_and_owned_output():
    m=mesh([[0,0],[4,0],[4,3],[0,3]])
    arrays=[x.copy() for x in (m.points,m.triangles,m.segments)]
    r=dm.improve(m,max_added_points=0)
    assert r.inserted_point is None and len(r.mesh.points)==4 and len(r.mesh.triangles)==2
    assert np.array_equal(r.mesh.points,m.points)
    for actual,expected in zip((m.points,m.triangles,m.segments),arrays):assert np.array_equal(actual,expected)
    r.mesh.points[0]=[99,99]
    assert np.array_equal(m.points,arrays[0])


@pytest.mark.parametrize('outer,holes',[
    ([[0,0],[8,0],[8,3],[4,3],[4,8],[0,8]],None),
    ([[0,0],[10,0],[10,10],[0,10]],[[[3,3],[3,7],[7,7],[7,3]]])])
def test_concavity_holes_and_fixed_internal_constraint(outer,holes):
    m=mesh(outer,holes)
    if holes is None:
        edges=np.sort(np.concatenate([m.triangles[:,[0,1]],m.triangles[:,[1,2]],m.triangles[:,[2,0]]]),axis=1)
        ids,counts=np.unique(edges,axis=0,return_counts=True)
        m.segments=np.vstack([m.segments,ids[counts==2][0]])
    r=dm.improve(m,max_added_points=1)
    assert area(r.mesh)==pytest.approx(area(m),rel=1e-12)
    assert np.array_equal(r.mesh.segments,m.segments)
    assert np.array_equal(r.segment_parents,np.arange(len(m.segments)))
    assert len(r.mesh.points)<=len(m.points)+1


def test_one_allowed_edge_has_oriented_parameter_and_segment_lineage():
    m,top=fan()
    r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert r.insertion_kind=='segment' and r.parent_segment==top
    assert r.inserted_point==len(m.points) and 0<r.segment_parameter<1
    a,b=m.points[m.segments[top]]
    assert np.allclose(r.mesh.points[-1],a+r.segment_parameter*(b-a),atol=1e-12)
    assert len(r.mesh.triangles)==len(m.triangles)+1
    children=r.mesh.segments[r.segment_parents==top]
    assert children.tolist()==[[int(m.segments[top,0]),len(m.points)],[len(m.points),int(m.segments[top,1])]]
    for i,s in enumerate(m.segments):
        if i!=top:assert r.mesh.segments[r.segment_parents==i].tolist()==[s.tolist()]
    assert r.after['penalty']<r.before['penalty']
    assert 1<=len(r.segment_trials)<=4
    assert r.segment_trials[0]==dict(parent_segment=top,segment_parameter=r.segment_parameter)
    assert all(x['parent_segment']==top and 0<x['segment_parameter']<1 for x in r.segment_trials)


def test_remote_unavoidable_triangle_does_not_veto_local_gain():
    m,top=fan()
    p=np.vstack([m.points,[[100,0],[101,0],[100,.00001]]])
    k=len(m.points)
    t=np.vstack([m.triangles,[[k,k+1,k+2]]])
    s=np.vstack([m.segments,[[k,k+1],[k+1,k+2],[k+2,k]]])
    combined=dm.Mesh(p,t,s,summary=None)
    r=dm.improve(combined,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert r.insertion_kind=='segment'
    assert r.after['minimum_q']==pytest.approx(r.before['minimum_q'])
    assert r.after['penalty']<r.before['penalty']
    assert any(set(x)=={k,k+1,k+2} for x in r.mesh.triangles)


def test_interior_permission_is_separate_from_segment_budget():
    m,top=fan(4,3)
    mixed=dm.improve(m,max_added_points=1,splittable_segments=[top])
    edge=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert mixed.insertion_kind=='interior' and edge.insertion_kind=='segment'
    assert mixed.after['penalty']<edge.after['penalty']<mixed.before['penalty']


def test_accepted_interior_then_boundary_split_preserves_both_points():
    m,top=fan(20,8)
    interior=dm.improve(m,max_added_points=1)
    assert interior.insertion_kind=='interior'
    edge=dm.improve(interior.mesh,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert edge.insertion_kind=='segment' and edge.parent_segment==top
    assert np.array_equal(edge.mesh.points[:-1],interior.mesh.points)
    assert interior.inserted_point in edge.mesh.triangles
    assert edge.inserted_point in edge.mesh.triangles
    # Retessellation receives every fixed point, including the previous interior.
    rebuilt=dm.mesh(dm.Domain(edge.mesh.points,edge.mesh.segments),
                    options=dm.MeshingOptions(refine=False,acute_protection='none'))
    assert np.array_equal(rebuilt.points,edge.mesh.points)
    assert set(np.unique(rebuilt.triangles))==set(range(len(edge.mesh.points)))
    assert area(rebuilt)==pytest.approx(area(m),rel=1e-12)


@pytest.mark.parametrize('offset,scale',[(np.zeros(2),1.),(np.array([386000.,6174000.]),1.),(np.array([-30.,60.]),.001),(np.zeros(2),1000.)])
def test_translation_scale_determinism_and_no_forbidden_split(offset,scale):
    m,top=fan();m.points=m.points*scale+offset
    r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    other=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert np.array_equal(r.mesh.points,other.mesh.points)
    assert np.array_equal(r.mesh.triangles,other.mesh.triangles)
    assert r.insertion_kind=='segment'
    assert np.array_equal(r.mesh.points[:len(m.points)],m.points)
    fixed=dm.improve(m,max_added_points=1,allow_interior=False)
    assert fixed.insertion_kind=='none' and len(fixed.mesh.points)==len(m.points)
    zero=dm.improve(m,max_added_points=0,splittable_segments=[top])
    assert len(zero.mesh.points)==len(m.points)


@pytest.mark.parametrize('budget',[-1,2,100,True,1.5])
def test_invalid_budgets(budget):
    m,_=fan()
    with pytest.raises(ValueError):dm.improve(m,max_added_points=budget)


def test_invalid_meshes_and_error_cleanup():
    m,_=fan()
    variants=[dm.Mesh(m.points,np.array([[0,1,99999]]),m.segments,None),
              dm.Mesh(m.points,m.triangles[:,::-1],m.segments,None),
              dm.Mesh(m.points,np.vstack([m.triangles,m.triangles[:1]]),m.segments,None),
              dm.Mesh(m.points,m.triangles,np.empty((0,2),int),None)]
    for _ in range(20):
        for invalid in variants:
            with pytest.raises(ValueError):dm.improve(invalid)
    with pytest.raises(ValueError):dm.improve(m,splittable_segments=[-1])
    assert dm.improve(m).mesh.triangles.size


def test_shape_flip_is_valid_and_explicitly_non_delaunay():
    p=np.array([[1.2084179604813374,1.913531611125363],[.4013367175598736,2.7131251356054227],
                [-1.7330339708859828,2.3773482795422525],[-.10840379431529797,-.8245773271285304]])
    m=mesh(p);r=dm.improve(m,max_added_points=0)
    assert r.flips==1 and r.after['minimum_q']>r.before['minimum_q']
    assert area(r.mesh)==pytest.approx(area(m)) and np.array_equal(r.mesh.points,m.points)
    def circle_test(triangle):
        other=next(i for i in range(4) if i not in triangle)
        q=p[triangle]-p[other]
        return np.linalg.det(np.c_[q,np.sum(q*q,axis=1)])
    assert circle_test(m.triangles[0])<=1e-10
    assert circle_test(r.mesh.triangles[0])>0


def test_moderate_fan_is_rewarded_without_severe_defects():
    m=mesh([[x,0] for x in (0,3,6,9,12)]+[[12,4],[0,4]])
    top=int(np.flatnonzero(np.all(m.points[m.segments,1]==4,axis=1))[0])
    r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert r.insertion_kind=='segment' and len(r.mesh.triangles)==len(m.triangles)+1
    assert r.quality_before.count_min_angle_lt_20==2
    assert r.mesh.summary.count_min_angle_lt_20==0
    assert r.before['poor_10_count']==r.after['poor_10_count']==0
    for value in (m,r.mesh):
        p=value.points[value.triangles];e=np.roll(p,-1,axis=1)-p
        a=abs(e[:,0,0]*e[:,1,1]-e[:,0,1]*e[:,1,0])/2
        q=4*np.sqrt(3)*a/(e*e).sum(axis=(1,2))
        # The severe-only criterion is blind to this moderate fan.
        assert np.maximum(0,.25-q).sum()==0
    assert r.before['penalty']-r.after['penalty']>.015


def test_segment_local_search_survives_remote_worst_triangle_selection():
    rng=np.random.default_rng(353)
    x=np.r_[0,np.sort(rng.uniform(.1,19.9,9)),20]
    m=mesh([[v,0] for v in x]+[[20,2],[0,2]])
    top=int(np.flatnonzero(np.all(m.points[m.segments,1]==2,axis=1))[0])
    for i in range(12):
        k=len(m.points);x=100+2*i
        m.points=np.vstack([m.points,[[x,0],[x+1,0],[x,.000001]]])
        m.triangles=np.vstack([m.triangles,[[k,k+1,k+2]]])
        m.segments=np.vstack([m.segments,[[k,k+1],[k+1,k+2],[k+2,k]]])
    r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
    assert r.insertion_kind=='segment'
    assert r.segment_parameter==pytest.approx(.2644311965128751)
    assert r.before['penalty']-r.after['penalty']>1.
    assert np.array_equal(r.mesh.points[:-1],m.points)
    assert r.after['poor_5_count']<=r.before['poor_5_count']
    assert r.after['poor_10_count']<=r.before['poor_10_count']


@pytest.mark.parametrize('offset,scale', [(np.zeros(2), 1.), (np.zeros(2), .001),
    (np.zeros(2), 1000.), (np.array([386000., 6174000.]), 1.)])
def test_thin_convex_hull_has_complete_coverage_and_fixed_constraints(offset, scale):
    import json
    from pathlib import Path
    data = json.loads((Path(__file__).parents[1]/'cases/thin_hull_domain.json').read_text())
    points = np.asarray(data['points'])*scale+offset
    segments = np.asarray(data['segments'], dtype=np.uint32)
    polygon = points[[0, 3, 4, 2, 1]]-points[0]
    expected = abs(np.sum(polygon[:, 0]*np.roll(polygon[:, 1], -1)
                          -polygon[:, 1]*np.roll(polygon[:, 0], -1)))/2
    raw = dm.mesh(dm.Domain(points, segments),
                  options=dm.MeshingOptions(refine=False, acute_protection='none'))
    assert len(raw.points) == 5 and len(raw.triangles) == 3
    assert np.array_equal(raw.points, points)
    assert set(map(tuple, np.sort(raw.segments, axis=1))) == set(map(tuple, np.sort(segments, axis=1)))
    assert area(raw) == pytest.approx(expected, rel=1e-10)
    assert area(dm.improve(raw, max_added_points=0).mesh) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize('case', ['fixed_interior_near_boundary', 'near_collinear_boundary'])
def test_nearly_collinear_boundary_with_fixed_interior_points(case):
    import json
    from pathlib import Path
    data = json.loads((Path(__file__).parents[1]/f'cases/{case}.json').read_text())
    points, segments = np.asarray(data['points']), np.asarray(data['segments'], dtype=np.uint32)
    raw = dm.mesh(dm.Domain(points, segments),
                  options=dm.MeshingOptions(refine=False, acute_protection='none'))
    assert np.array_equal(raw.points, points)
    assert set(map(tuple, np.sort(raw.segments, axis=1))) == set(map(tuple, np.sort(segments, axis=1)))
    boundary_points = len(np.unique(segments))
    assert len(raw.triangles) == 2*len(points)-boundary_points-2
    assert set(np.unique(raw.triangles)) == set(range(len(points)))
    assert area(dm.improve(raw, max_added_points=0).mesh) == pytest.approx(area(raw), rel=1e-10)
