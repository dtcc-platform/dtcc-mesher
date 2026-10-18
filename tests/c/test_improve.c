#include "dtcc_mesher/dtcc_mesher_improve.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define REQUIRE(x) do {if(!(x)){fprintf(stderr,"failed line %d: %s\n",__LINE__,#x);return 1;}}while(0)
int main(void)
{
    dtcc_mesher_point p[]={{0,0},{4,0},{4,3},{0,3}};
    uint32_t triangles[]={0,1,2,0,2,3};
    dtcc_mesher_segment segments[]={{0,1},{1,2},{2,3},{3,0}};
    dtcc_mesher_mesh mesh={0};dtcc_mesher_improvement result;dtcc_mesher_error error;
    mesh.points=p;mesh.num_points=4;mesh.triangles=triangles;mesh.num_triangles=2;mesh.segments=segments;mesh.num_segments=4;
    dtcc_mesher_improvement_init(&result);
    for(int i=0;i<100;i++) {
        REQUIRE(dtcc_mesher_improve(&mesh,0,1,NULL,0,&result,&error)==DTCC_MESHER_STATUS_OK);
        REQUIRE(result.mesh.num_points==4 && result.mesh.num_triangles==2 && result.insertion_kind==DTCC_MESHER_INSERT_NONE);
        REQUIRE(result.num_segment_trials==0);
        REQUIRE(result.mesh.points!=mesh.points && !memcmp(result.mesh.points,p,sizeof(p)));
        for(size_t j=0;j<4;j++) REQUIRE(result.segment_parents[j]==j);
        REQUIRE(dtcc_mesher_improve(&result.mesh,0,1,NULL,0,&result,&error)==DTCC_MESHER_STATUS_INVALID_ARGUMENT);
        REQUIRE(result.mesh.num_points==4);
        dtcc_mesher_improvement_free(&result);dtcc_mesher_improvement_free(&result);
        REQUIRE(!result.mesh.points && !result.segment_parents);
        REQUIRE(dtcc_mesher_improve(&mesh,2,1,NULL,0,&result,&error)==DTCC_MESHER_STATUS_INVALID_ARGUMENT);
        uint32_t saved=triangles[0];triangles[0]=100;
        REQUIRE(dtcc_mesher_improve(&mesh,0,1,NULL,0,&result,&error)==DTCC_MESHER_STATUS_GEOMETRY);
        REQUIRE(!result.mesh.points && !result.segment_parents);triangles[0]=saved;
    }
    dtcc_mesher_improvement_free(NULL);
    dtcc_mesher_point fan[23];dtcc_mesher_segment edges[23];dtcc_mesher_domain domain={0};dtcc_mesher_options options;
    dtcc_mesher_mesh generated={0};
    for(int i=0;i<21;i++){fan[i].x=i;fan[i].y=0;}
    fan[21]=(dtcc_mesher_point){20,8};fan[22]=(dtcc_mesher_point){0,8};
    for(int i=0;i<23;i++) edges[i]=(dtcc_mesher_segment){(uint32_t)i,(uint32_t)((i+1)%23)};
    domain.points=fan;domain.num_points=23;domain.segments=edges;domain.num_segments=23;
    dtcc_mesher_options_init(&options);options.enable_refinement=0;options.enable_acute_protection=0;
    REQUIRE(dtcc_mesher_generate(&domain,&options,&generated,&error)==DTCC_MESHER_STATUS_OK);
    size_t top=SIZE_MAX;
    for(size_t i=0;i<generated.num_segments;i++) if(generated.points[generated.segments[i].a].y==8 && generated.points[generated.segments[i].b].y==8) top=i;
    REQUIRE(top!=SIZE_MAX);
    REQUIRE(dtcc_mesher_improve(&generated,1,0,&top,1,&result,&error)==DTCC_MESHER_STATUS_OK);
    REQUIRE(result.insertion_kind==DTCC_MESHER_INSERT_SEGMENT && result.parent_segment==top);
    REQUIRE(result.mesh.num_points==generated.num_points+1 && result.mesh.num_triangles==generated.num_triangles+1);
    REQUIRE(fabs(result.segment_parameter-.5)<1e-12 && result.after.penalty<result.before.penalty);
    REQUIRE(result.num_segment_trials>=1 && result.num_segment_trials<=DTCC_MESHER_SEGMENT_TRIAL_LIMIT);
    REQUIRE(result.segment_trials[0].parent_segment==top);
    REQUIRE(result.segment_trials[0].segment_parameter==result.segment_parameter);
    for(size_t i=0;i<result.num_segment_trials;i++) {
        REQUIRE(result.segment_trials[i].parent_segment==top);
        REQUIRE(result.segment_trials[i].segment_parameter>0 && result.segment_trials[i].segment_parameter<1);
    }
    dtcc_mesher_improvement_free(&result);dtcc_mesher_mesh_free(&generated);
    return 0;
}
