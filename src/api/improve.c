#include "dtcc_mesher/dtcc_mesher_improve.h"
#include "dtcc_mesher_api_internal.h"
#include "cdt.h"
#include "validate.h"
#include "predicates.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CANDIDATE_LIMIT 128
#define POOR_TRIANGLE_LIMIT 12
#define SEGMENT_TRIANGLE_LIMIT 6
#define FLIP_LIMIT 4096

static double orientation(const double *a, const double *b, const double *c)
{ return orient2d((REAL *)a, (REAL *)b, (REAL *)c); }

static double metric(const TMMesh *m, const TMTriangle *t, double *area, double *angle)
{
    const double *p[3]; double sum=0., minimum=180.; int k;
    for(k=0;k<3;k++) p[k]=m->points[t->v[k]].xy;
    *area=orientation(p[0],p[1],p[2])/2.;
    for(k=0;k<3;k++) {
        double ux=p[(k+1)%3][0]-p[k][0], uy=p[(k+1)%3][1]-p[k][1];
        double vx=p[(k+2)%3][0]-p[k][0], vy=p[(k+2)%3][1]-p[k][1];
        double a=atan2(fabs(ux*vy-uy*vx),ux*vx+uy*vy)*57.29577951308232;
        sum+=ux*ux+uy*uy; if(a<minimum) minimum=a;
    }
    *angle=minimum;
    return sum>0 ? 4.*sqrt(3.)*(*area)/sum : 0.;
}

static dtcc_mesher_shape_score score(const TMMesh *m)
{
    dtcc_mesher_shape_score s={0}; size_t i; s.minimum_q=1.;
    for(i=0;i<m->triangle_count;i++) {
        double area, angle, q=metric(m,&m->triangles[i],&area,&angle);
        s.penalty+=(1.-q)*(1.-q); s.minimum_q=fmin(s.minimum_q,q);
        if(angle<5.-1e-7) {s.poor_5_count++;s.poor_5_area+=area;}
        if(angle<10.-1e-7) {s.poor_10_count++;s.poor_10_area+=area;}
    }
    return s;
}

static int guarded(dtcc_mesher_shape_score a, dtcc_mesher_shape_score b)
{
    return b.minimum_q>=a.minimum_q-1e-10 && b.poor_5_count<=a.poor_5_count &&
        b.poor_5_area<=a.poor_5_area+1e-10*fmax(1.,a.poor_5_area) &&
        b.poor_10_count<=a.poor_10_count &&
        b.poor_10_area<=a.poor_10_area+1e-10*fmax(1.,a.poor_10_area);
}

static TMStatus clone(const TMMesh *a,TMMesh *b)
{
    memset(b,0,sizeof(*b));
    b->point_count=b->point_capacity=a->point_count;
    b->triangle_count=b->triangle_capacity=a->triangle_count;
    b->segment_count=b->segment_capacity=a->segment_count;
    b->points=malloc(b->point_count*sizeof(*b->points));
    b->triangles=malloc(b->triangle_count*sizeof(*b->triangles));
    b->segments=malloc(b->segment_count*sizeof(*b->segments));
    if(!b->points || !b->triangles || (b->segment_count && !b->segments)) {tm_free_mesh(b);return TM_ERR_ALLOC;}
    memcpy(b->points,a->points,b->point_count*sizeof(*b->points));
    memcpy(b->triangles,a->triangles,b->triangle_count*sizeof(*b->triangles));
    if(b->segment_count) memcpy(b->segments,a->segments,b->segment_count*sizeof(*b->segments));
    return TM_OK;
}

/* Positive orientation and a conforming planar embedding establish disjoint
 * triangle interiors. Check crossings, containment and T junctions on import.
 * Subsequent changes use the native legal topology primitives only.
 */
static TMStatus validate_embedding(const TMMesh *m)
{
    size_t i,j; int e,f;
    for(i=0;i<m->triangle_count;i++) {
        const TMTriangle *t=&m->triangles[i];
        for(j=0;j<m->point_count;j++) {
            double a[3]; int k;
            if(tm_triangle_contains_vertex(t,(int)j)) continue;
            for(k=0;k<3;k++) a[k]=orientation(m->points[t->v[k]].xy,m->points[t->v[(k+1)%3]].xy,m->points[j].xy);
            if(a[0]>0 && a[1]>0 && a[2]>0) return TM_ERR_INVALID_MESH;
            for(k=0;k<3;k++) if(a[k]==0.) {
                const double *x=m->points[t->v[k]].xy,*y=m->points[t->v[(k+1)%3]].xy,*p=m->points[j].xy;
                if((p[0]-x[0])*(p[0]-y[0])+(p[1]-x[1])*(p[1]-y[1])<0.) return TM_ERR_INVALID_MESH;
            }
        }
        for(j=i+1;j<m->triangle_count;j++) for(e=0;e<3;e++) for(f=0;f<3;f++) {
            const double *a=m->points[t->v[e]].xy,*b=m->points[t->v[(e+1)%3]].xy;
            const double *c=m->points[m->triangles[j].v[f]].xy,*d=m->points[m->triangles[j].v[(f+1)%3]].xy;
            double x=orientation(a,b,c),y=orientation(a,b,d),u=orientation(c,d,a),v=orientation(c,d,b);
            if(((x<0 && y>0)||(x>0 && y<0)) && ((u<0 && v>0)||(u>0 && v<0))) return TM_ERR_INVALID_MESH;
        }
    }
    return TM_OK;
}

static TMStatus import_mesh(const dtcc_mesher_mesh *a,TMMesh *b,double offset[2])
{
    size_t i,j; TMStatus status; TMValidationReport report;
    memset(b,0,sizeof(*b));
    if(!a || !a->points || !a->triangles || a->num_points<3 || !a->num_triangles ||
       a->num_points>4096 || a->num_triangles>4096 || a->num_segments>3*a->num_triangles ||
       (a->num_segments && !a->segments)) return TM_ERR_INVALID_MESH;
    b->point_count=b->point_capacity=a->num_points;
    b->triangle_count=b->triangle_capacity=a->num_triangles;
    b->segment_count=b->segment_capacity=a->num_segments;
    b->points=calloc(b->point_count,sizeof(*b->points));
    b->triangles=calloc(b->triangle_count,sizeof(*b->triangles));
    b->segments=calloc(b->segment_count,sizeof(*b->segments));
    if(!b->points || !b->triangles || (b->segment_count && !b->segments)) return TM_ERR_ALLOC;
    offset[0]=a->points[0].x;offset[1]=a->points[0].y;
    for(i=0;i<b->point_count;i++) {
        TMPoint *p=&b->points[i];
        p->xy[0]=a->points[i].x-offset[0];p->xy[1]=a->points[i].y-offset[1];
        if(!isfinite(p->xy[0]) || !isfinite(p->xy[1])) return TM_ERR_INVALID_MESH;
        p->original_index=(int)i;p->kind=TM_VERTEX_INPUT;p->incident_triangle=-1;p->protection_apex=-1;
        for(j=0;j<i;j++) if(p->xy[0]==b->points[j].xy[0] && p->xy[1]==b->points[j].xy[1]) return TM_ERR_DUPLICATE;
    }
    for(i=0;i<b->triangle_count;i++) {
        int k;
        for(k=0;k<3;k++) {
            uint32_t v=a->triangles[3*i+k];if(v>=b->point_count) return TM_ERR_INVALID_MESH;
            b->triangles[i].v[k]=(int)v;
        }
        if(orientation(b->points[b->triangles[i].v[0]].xy,b->points[b->triangles[i].v[1]].xy,
                       b->points[b->triangles[i].v[2]].xy)<=0.) return TM_ERR_INVALID_MESH;
    }
    for(i=0;i<b->segment_count;i++) {
        TMSegment *s=&b->segments[i];
        if(a->segments[i].a>=b->point_count || a->segments[i].b>=b->point_count || a->segments[i].a==a->segments[i].b) return TM_ERR_INVALID_MESH;
        s->v[0]=(int)a->segments[i].a;s->v[1]=(int)a->segments[i].b;s->original_index=(int)i;s->live=1;s->protected_apex=-1;
        for(j=0;j<i;j++) if((s->v[0]==b->segments[j].v[0] && s->v[1]==b->segments[j].v[1]) ||
                            (s->v[0]==b->segments[j].v[1] && s->v[1]==b->segments[j].v[0])) return TM_ERR_INVALID_MESH;
    }
    status=tm_rebuild_topology(b);if(status!=TM_OK) return status;
    for(i=0;i<b->triangle_count;i++) {
        int e;
        for(e=0;e<3;e++) {
            int x,y,found=0;tm_triangle_edge_vertices(&b->triangles[i],e,&x,&y);
            for(j=0;j<b->segment_count;j++) if(tm_find_edge_in_triangle(&b->triangles[i],b->segments[j].v[0],b->segments[j].v[1])==e) found=1;
            b->triangles[i].constrained[e]=(unsigned char)found;
            if(b->triangles[i].nbr[e]<0 && !found) return TM_ERR_INVALID_MESH;
            if(b->triangles[i].nbr[e]>=0) {
                int n=b->triangles[i].nbr[e],ne=tm_find_edge_in_triangle(&b->triangles[n],x,y);
                if(orientation(b->points[x].xy,b->points[y].xy,b->points[b->triangles[n].v[ne]].xy)>=0.) return TM_ERR_INVALID_MESH;
            }
        }
    }
    for(i=0;i<b->point_count;i++) if(b->points[i].incident_triangle<0) return TM_ERR_INVALID_MESH;
    status=tm_validate_mesh(b,0,&report);if(status!=TM_OK) return status;
    return validate_embedding(b);
}

static TMStatus flips(TMMesh *m,size_t *count)
{
    size_t limit=8*m->triangle_count+1,round; if(limit>FLIP_LIMIT) limit=FLIP_LIMIT;
    for(round=0;round<limit;round++) {
        size_t i;int changed=0;
        for(i=0;i<m->triangle_count && !changed;i++) {
            int e;
            for(e=0;e<3;e++) {
                int n=m->triangles[i].nbr[e],a,b,ne; TMTriangle old[2], proposed[2];
                TMMesh pair; double area,angle,q0,q1,r0,r1;
                dtcc_mesher_shape_score before,after;
                if(n<(int)i || !tm_edge_is_flippable(m,(int)i,e)) continue;
                old[0]=m->triangles[i];old[1]=m->triangles[n];
                tm_triangle_edge_vertices(&old[0],e,&a,&b);ne=tm_find_edge_in_triangle(&old[1],a,b);
                proposed[0]=old[0];proposed[1]=old[1];
                proposed[0].v[0]=old[0].v[e];proposed[0].v[1]=old[1].v[ne];proposed[0].v[2]=b;
                proposed[1].v[0]=old[1].v[ne];proposed[1].v[1]=old[0].v[e];proposed[1].v[2]=a;
                /* q is orientation independent here; topology is delegated to tm_flip_edge. */
                q0=metric(m,&old[0],&area,&angle);q1=metric(m,&old[1],&area,&angle);
                r0=fabs(metric(m,&proposed[0],&area,&angle));r1=fabs(metric(m,&proposed[1],&area,&angle));
                if(fmin(r0,r1)<=fmin(q0,q1)+1e-10) continue;
                pair=*m;pair.triangle_count=2;pair.triangles=old;before=score(&pair);
                /* Use the native flip to evaluate its exact oriented result. */
                {TMStatus s=tm_flip_edge(m,(int)i,e);if(s!=TM_OK) return s;}
                proposed[0]=m->triangles[i];proposed[1]=m->triangles[n];pair.triangles=proposed;after=score(&pair);
                if(after.penalty>before.penalty+1e-12 || !guarded(before,after)) {
                    int undo=tm_find_edge_in_triangle(&m->triangles[i],old[0].v[e],old[1].v[ne]);
                    TMStatus s=tm_flip_edge(m,(int)i,undo);if(s!=TM_OK) return s;
                    continue;
                }
                (*count)++;changed=1;break;
            }
        }
        if(!changed) break;
    }
    return TM_OK;
}

static int same_triangle(const TMTriangle *a,const TMTriangle *b)
{ return tm_triangle_contains_vertex(a,b->v[0]) && tm_triangle_contains_vertex(a,b->v[1]) && tm_triangle_contains_vertex(a,b->v[2]); }

static double changed_minimum(const TMMesh *a,const TMMesh *b)
{
    size_t i,j;double q=1.;
    for(i=0;i<a->triangle_count;i++) {
        int found=0;double area,angle;
        for(j=0;j<b->triangle_count;j++) if(same_triangle(&a->triangles[i],&b->triangles[j])) {found=1;break;}
        if(!found) q=fmin(q,metric(a,&a->triangles[i],&area,&angle));
    }
    return q;
}

typedef struct {double p[2],t;int segment;} Candidate;
static void propose(Candidate *c,size_t *count,const TMMesh *m,double x,double y,int segment,double t)
{
    size_t i;double scale=fmax(1.,fmax(fabs(x),fabs(y))),tol=64*DBL_EPSILON*scale;
    if(!isfinite(x) || !isfinite(y) || *count>=CANDIDATE_LIMIT) return;
    for(i=0;i<m->point_count;i++) if(hypot(x-m->points[i].xy[0],y-m->points[i].xy[1])<=tol) return;
    for(i=0;i<*count;i++) if(hypot(x-c[i].p[0],y-c[i].p[1])<=tol) return;
    c[*count].p[0]=x;c[*count].p[1]=y;c[*count].segment=segment;c[*count].t=t;(*count)++;
}

static size_t candidates(const TMMesh *m,const size_t *allowed,size_t n,int interior,Candidate *out)
{
    size_t count=0,i,j,k,selected[POOR_TRIANGLE_LIMIT],ns=0;double qualities[POOR_TRIANGLE_LIMIT];
    for(i=0;i<m->triangle_count;i++) {
        double area,angle,q=metric(m,&m->triangles[i],&area,&angle);
        for(j=0;j<ns;j++) if(q<qualities[j]) break;
        if(j>=POOR_TRIANGLE_LIMIT) continue;
        if(ns<POOR_TRIANGLE_LIMIT) ns++;
        for(k=ns-1;k>j;k--) {selected[k]=selected[k-1];qualities[k]=qualities[k-1];}
        selected[j]=i;qualities[j]=q;
    }
    /* Segment proposals use one midpoint and projections of poor triangle
     * centroids/opposite vertices. They are feature-informed, never a grid.
     */
    for(i=0;i<n;i++) {
        const TMSegment *s=&m->segments[allowed[i]];
        const double *a=m->points[s->v[0]].xy,*b=m->points[s->v[1]].xy;
        double dx=b[0]-a[0],dy=b[1]-a[1],l2=dx*dx+dy*dy;
        propose(out,&count,m,a[0]+.5*dx,a[1]+.5*dy,(int)allowed[i],.5);
        /* A remote acute component must not consume this segment's search.
         * Select a bounded endpoint star by q, with stable triangle-index ties.
         * These are local shape proposals, not opposite-boundary samples.
         */
        {
            size_t local[SEGMENT_TRIANGLE_LIMIT],nl=0,h;
            double local_q[SEGMENT_TRIANGLE_LIMIT];
            for(j=0;j<m->triangle_count;j++) {
                const TMTriangle *tri=&m->triangles[j];
                double area,angle,q;
                if(!tm_triangle_contains_vertex(tri,s->v[0]) &&
                   !tm_triangle_contains_vertex(tri,s->v[1])) continue;
                q=metric(m,tri,&area,&angle);
                for(k=0;k<nl;k++) if(q<local_q[k]) break;
                if(k>=SEGMENT_TRIANGLE_LIMIT) continue;
                if(nl<SEGMENT_TRIANGLE_LIMIT) nl++;
                for(h=nl-1;h>k;h--) {local[h]=local[h-1];local_q[h]=local_q[h-1];}
                local[k]=j;local_q[k]=q;
            }
            for(j=0;j<nl;j++) {
                const TMTriangle *tri=&m->triangles[local[j]];
                double x=0.,y=0.,t;
                for(k=0;k<3;k++) {x+=m->points[tri->v[k]].xy[0]/3.;y+=m->points[tri->v[k]].xy[1]/3.;}
                t=((x-a[0])*dx+(y-a[1])*dy)/l2;
                if(t>.02 && t<.98) propose(out,&count,m,a[0]+t*dx,a[1]+t*dy,(int)allowed[i],t);
                if(tm_triangle_contains_vertex(tri,s->v[0]) && tm_triangle_contains_vertex(tri,s->v[1])) {
                    for(k=0;k<3;k++) if(tri->v[k]!=s->v[0] && tri->v[k]!=s->v[1]) {
                        const double *p=m->points[tri->v[k]].xy;
                        t=((p[0]-a[0])*dx+(p[1]-a[1])*dy)/l2;
                        if(t>.02 && t<.98) propose(out,&count,m,a[0]+t*dx,a[1]+t*dy,(int)allowed[i],t);
                    }
                }
            }
        }
        for(j=0;j<ns;j++) {
            const TMTriangle *tri=&m->triangles[selected[j]];
            double x=0.,y=0.,t;
            for(k=0;k<3;k++) {x+=m->points[tri->v[k]].xy[0]/3.;y+=m->points[tri->v[k]].xy[1]/3.;}
            t=((x-a[0])*dx+(y-a[1])*dy)/l2;
            if(t>.02 && t<.98) propose(out,&count,m,a[0]+t*dx,a[1]+t*dy,(int)allowed[i],t);
            if(tm_triangle_contains_vertex(tri,s->v[0]) && tm_triangle_contains_vertex(tri,s->v[1])) {
                for(k=0;k<3;k++) if(tri->v[k]!=s->v[0] && tri->v[k]!=s->v[1]) {
                    const double *p=m->points[tri->v[k]].xy;t=((p[0]-a[0])*dx+(p[1]-a[1])*dy)/l2;
                    if(t>.02 && t<.98) propose(out,&count,m,a[0]+t*dx,a[1]+t*dy,(int)allowed[i],t);
                }
            }
        }
    }
    for(i=0;interior && i<ns;i++) {
        const TMTriangle *t=&m->triangles[selected[i]];double x=0.,y=0.;
        for(k=0;k<3;k++) {x+=m->points[t->v[k]].xy[0]/3.;y+=m->points[t->v[k]].xy[1]/3.;}
        propose(out,&count,m,x,y,-1,0.);
        for(k=0;k<3;k++) {
            const double *a=m->points[t->v[k]].xy,*b=m->points[t->v[(k+1)%3]].xy,*d=m->points[t->v[(k+2)%3]].xy;
            double dx=b[0]-a[0],dy=b[1]-a[1],sign=dx*(d[1]-a[1])-dy*(d[0]-a[0]);
            double factor=sign>0 ? .866025403784 : -.866025403784;
            propose(out,&count,m,(a[0]+b[0])/2.-factor*dy,(a[1]+b[1])/2.+factor*dx,-1,0.);
        }
    }
    return count;
}

void dtcc_mesher_improvement_init(dtcc_mesher_improvement *r)
{ if(r) {memset(r,0,sizeof(*r));r->parent_segment=SIZE_MAX;r->inserted_point=SIZE_MAX;} }
void dtcc_mesher_improvement_free(dtcc_mesher_improvement *r)
{ if(r) {dtcc_mesher_mesh_free(&r->mesh);free(r->segment_parents);dtcc_mesher_improvement_init(r);} }

dtcc_mesher_status dtcc_mesher_improve(const dtcc_mesher_mesh *input,unsigned int budget,int interior,
    const size_t *allowed,size_t n,dtcc_mesher_improvement *r,dtcc_mesher_error *error)
{
    TMMesh base={0},best={0};TMStatus s=TM_OK;TMValidationReport validation;
    Candidate proposals[CANDIDATE_LIMIT];size_t count=0,i,j,best_flips=0,base_flips=0,best_index=SIZE_MAX;
    double segment_benefits[DTCC_MESHER_SEGMENT_TRIAL_LIMIT]={0};
    double offset[2],best_benefit=0.;const char *stage="input";dtcc_mesher_shape_score baseline;dtcc_mesher_status public_status;
    if(!r || (input && input==&r->mesh) || budget>1 || (interior!=0 && interior!=1) || (n && !allowed)) {
        dtcc_mesher_api_set_error(error,DTCC_MESHER_STATUS_INVALID_ARGUMENT,"improvement budget must be zero or one; output and allowed-segment arrays must be valid");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }
    dtcc_mesher_improvement_init(r);dtcc_mesher_api_clear_error(error);tm_initialize();
    s=import_mesh(input,&base,offset);if(s!=TM_OK) goto fail;
    if(n>base.segment_count) {s=TM_ERR_INVALID_PSLG;goto fail;}
    for(i=0;i<n;i++) {
        if(allowed[i]>=base.segment_count) {s=TM_ERR_INVALID_PSLG;goto fail;}
        for(j=0;j<i;j++) if(allowed[i]==allowed[j]) {s=TM_ERR_INVALID_PSLG;goto fail;}
    }
    stage="zero_point_optimization";r->before=score(&base);
    s=flips(&base,&base_flips);if(s!=TM_OK) goto fail;
    baseline=score(&base);s=clone(&base,&best);if(s!=TM_OK) goto fail;
    stage="one_point_trial";
    if(budget) count=candidates(&base,allowed,n,interior,proposals);
    for(i=0;i<count;i++) {
        TMMesh trial={0};int point=-1,attempted=0;size_t nf=0;TMLocation location;dtcc_mesher_shape_score after;
        double benefit,cost,total_area0=0.,total_area1=0.;
        s=clone(&base,&trial);if(s!=TM_OK) goto fail;
        r->candidate_evaluations++;
        if(proposals[i].segment>=0) {
            int found=0;const TMSegment *seg=&base.segments[proposals[i].segment];
            for(j=0;j<trial.triangle_count;j++) {
                int e=tm_find_edge_in_triangle(&trial.triangles[j],seg->v[0],seg->v[1]);
                if(e>=0) {
                    found=1;attempted=1;s=tm_insert_point_on_edge(&trial,(int)j,e,proposals[i].p,TM_VERTEX_SEGMENT_SPLIT,&point);break;
                }
            }
            if(!found) s=TM_ERR_INVALID_MESH;
            if(s==TM_OK) s=tm_split_segment_registry(&trial,(size_t)proposals[i].segment,point);
        } else {
            s=tm_locate_point(&trial,proposals[i].p,0,&location);
            /* Location can return a boundary triangle for an exterior point.
             * Confirm strict containment; also cover disconnected components.
             */
            location.triangle=-1;
            for(j=0;j<trial.triangle_count;j++) {
                const TMTriangle *t=&trial.triangles[j];int k,inside=1;
                for(k=0;k<3;k++) if(orientation(trial.points[t->v[k]].xy,
                    trial.points[t->v[(k+1)%3]].xy,proposals[i].p)<=0.) inside=0;
                if(inside) {location.triangle=(int)j;break;}
            }
            if(location.triangle>=0) {
                attempted=1;s=tm_insert_point_in_triangle(&trial,location.triangle,proposals[i].p,TM_VERTEX_TRIANGLE_SPLIT,&point);
            } else s=TM_ERR_INVALID_PSLG;
        }
        if(s==TM_ERR_ALLOC || (attempted && s!=TM_OK)) {tm_free_mesh(&trial);goto fail;}
        if(s!=TM_OK || trial.point_count!=base.point_count+1) {tm_free_mesh(&trial);continue;}
        s=flips(&trial,&nf);if(s!=TM_OK) {tm_free_mesh(&trial);goto fail;}
        s=tm_validate_mesh(&trial,0,&validation);if(s!=TM_OK) {tm_free_mesh(&trial);goto fail;}
        for(j=0;j<base.triangle_count;j++) {double area,angle;metric(&base,&base.triangles[j],&area,&angle);total_area0+=area;}
        for(j=0;j<trial.triangle_count;j++) {double area,angle;metric(&trial,&trial.triangles[j],&area,&angle);total_area1+=area;}
        if(fabs(total_area0-total_area1)>1e-10*fmax(1.,total_area0)) {tm_free_mesh(&trial);s=TM_ERR_INVALID_MESH;goto fail;}
        after=score(&trial);cost=(double)(trial.triangle_count-base.triangle_count);
        benefit=(baseline.penalty-after.penalty)/cost;
        if(guarded(baseline,after) && changed_minimum(&trial,&base)>=changed_minimum(&base,&trial)-1e-10 && benefit>.015) {
            if(proposals[i].segment>=0) {
                size_t rank,h;
                for(rank=0;rank<r->num_segment_trials;rank++)
                    if(benefit>segment_benefits[rank]+1e-12) break;
                if(rank<DTCC_MESHER_SEGMENT_TRIAL_LIMIT) {
                    if(r->num_segment_trials<DTCC_MESHER_SEGMENT_TRIAL_LIMIT) r->num_segment_trials++;
                    for(h=r->num_segment_trials-1;h>rank;h--) {
                        r->segment_trials[h]=r->segment_trials[h-1];segment_benefits[h]=segment_benefits[h-1];
                    }
                    r->segment_trials[rank].parent_segment=(size_t)proposals[i].segment;
                    r->segment_trials[rank].segment_parameter=proposals[i].t;segment_benefits[rank]=benefit;
                }
            }
            if(benefit>best_benefit+1e-12) {
                tm_free_mesh(&best);best=trial;memset(&trial,0,sizeof(trial));
                best_index=i;best_benefit=benefit;best_flips=nf;
            }
        }
        tm_free_mesh(&trial);
    }
    stage="output_validation";s=tm_validate_mesh(&best,0,&validation);if(s!=TM_OK) goto fail;
    r->after=score(&best);r->flips=base_flips+best_flips;
    public_status=dtcc_mesher_api_copy_internal_mesh(&best,&r->mesh,error);
    if(public_status!=DTCC_MESHER_STATUS_OK) {s=TM_ERR_ALLOC;goto fail;}
    for(i=0;i<r->mesh.num_points;i++) {
        if(i<input->num_points) r->mesh.points[i]=input->points[i];
        else {r->mesh.points[i].x+=offset[0];r->mesh.points[i].y+=offset[1];}
    }
    if(best_index!=SIZE_MAX) {
        const Candidate *c=&proposals[best_index];r->inserted_point=input->num_points;
        r->insertion_kind=c->segment<0 ? DTCC_MESHER_INSERT_INTERIOR : DTCC_MESHER_INSERT_SEGMENT;
        if(c->segment>=0) {r->parent_segment=(size_t)c->segment;r->segment_parameter=c->t;}
    }
    r->segment_parents=calloc(r->mesh.num_segments,sizeof(*r->segment_parents));
    if(r->mesh.num_segments && !r->segment_parents) {s=TM_ERR_ALLOC;goto fail;}
    for(i=0,j=0;i<best.segment_count;i++) if(best.segments[i].live) r->segment_parents[j++]=(size_t)best.segments[i].original_index;
    public_status=dtcc_mesher_analyze_mesh(input,&r->quality_before,error);
    if(public_status==DTCC_MESHER_STATUS_OK) public_status=dtcc_mesher_analyze_mesh(&r->mesh,&r->quality_after,error);
    if(public_status!=DTCC_MESHER_STATUS_OK) {s=TM_ERR_INTERNAL;goto fail;}
    snprintf(r->termination,sizeof(r->termination),"%s",best_index!=SIZE_MAX ? "accepted_one_point" : base_flips ? "zero_point_improved" : budget ? "no_beneficial_trial" : "zero_point_complete");
    tm_free_mesh(&base);tm_free_mesh(&best);return DTCC_MESHER_STATUS_OK;
fail:
    tm_free_mesh(&base);tm_free_mesh(&best);dtcc_mesher_improvement_free(r);
    public_status=(s==TM_ERR_ALLOC || !strcmp(stage,"input")) ? dtcc_mesher_api_map_status(s) : DTCC_MESHER_STATUS_INTERNAL;
    dtcc_mesher_api_set_error(error,public_status,"mesh improvement (%s): %s",stage,tm_internal_status_string(s));return public_status;
}
