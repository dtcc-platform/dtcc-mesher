#ifndef DTCC_MESHER_MESH_H
#define DTCC_MESHER_MESH_H

#include <stddef.h>

typedef enum {
    TM_VERTEX_INPUT = 0,
    TM_VERTEX_SEGMENT_SPLIT,
    TM_VERTEX_TRIANGLE_SPLIT,
    TM_VERTEX_SUPER
} TMVertexKind;

typedef enum {
    TM_PROTECTION_SIDE_NONE = 0,
    TM_PROTECTION_SIDE_LEFT,
    TM_PROTECTION_SIDE_RIGHT
} TMProtectionSide;

typedef struct {
    double xy[2];
    int original_index;
    TMVertexKind kind;
    int incident_triangle;
    int protection_apex;
    unsigned char protection_side;
    unsigned int protection_level;
} TMPoint;

/* Edge i is opposite vertex v[i]. */
typedef struct {
    int v[3];
    int nbr[3];
    unsigned char constrained[3];
    int dead;
} TMTriangle;

typedef struct {
    int v[2];
    int original_index;
    int live;
    int is_protected;
    int protected_apex;
} TMSegment;

typedef struct {
    TMPoint *points;
    size_t point_count;
    size_t point_capacity;
    TMTriangle *triangles;
    size_t triangle_count;
    size_t triangle_capacity;
    TMSegment *segments;
    size_t segment_count;
    size_t segment_capacity;
} TMMesh;

typedef struct {
    int triangle;
    int edge;
    int on_edge;
} TMLocation;

typedef enum {
    TM_OK = 0,
    TM_ERR_ALLOC,
    TM_ERR_IO,
    TM_ERR_PARSE,
    TM_ERR_DUPLICATE,
    TM_ERR_TOO_FEW_POINTS,
    TM_ERR_COLLINEAR,
    TM_ERR_INVALID_PSLG,
    TM_ERR_INVALID_MESH,
    TM_ERR_INTERNAL
} TMStatus;

void tm_initialize(void);
const char *tm_internal_status_string(TMStatus status);
TMStatus tm_read_points_file(const char *path, TMPoint **out_points, size_t *out_count);
TMStatus tm_build_mesh(const TMPoint *input_points, size_t input_count, TMMesh *out_mesh);
TMStatus tm_rebuild_topology(TMMesh *mesh);
void tm_triangle_edge_vertices(const TMTriangle *triangle, int edge, int *out_a, int *out_b);
int tm_find_edge_in_triangle(const TMTriangle *triangle, int a, int b);
int tm_triangle_contains_vertex(const TMTriangle *triangle, int vertex_index);
int tm_edge_is_locally_delaunay(const TMMesh *mesh, int triangle_index, int edge);
int tm_triangle_is_quality_exempt(const TMMesh *mesh, int triangle_index);
size_t tm_count_protected_corners(const TMMesh *mesh);
TMStatus tm_locate_point(const TMMesh *mesh, const double point[2], int start_triangle, TMLocation *out_location);
TMStatus tm_flip_edge(TMMesh *mesh, int triangle_index, int edge);
TMStatus tm_insert_point_in_triangle(
    TMMesh *mesh,
    int triangle_index,
    const double point[2],
    TMVertexKind kind,
    int *out_point_index
);
TMStatus tm_insert_point_on_edge(
    TMMesh *mesh,
    int triangle_index,
    int edge,
    const double point[2],
    TMVertexKind kind,
    int *out_point_index
);
void tm_free_points(TMPoint *points);
void tm_free_mesh(TMMesh *mesh);

#endif
