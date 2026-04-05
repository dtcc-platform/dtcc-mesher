#include "mesh.h"

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_pslg.h"
#include "predicates.h"

typedef struct {
    int v[3];
    int dead;
} TMWorkingTriangle;

typedef struct {
    int a;
    int b;
    int active;
} TMEdge;

typedef struct {
    int a;
    int b;
} TMEdgePair;

typedef struct {
    int triangle_index;
    int edge;
} TMFlipCandidate;

static int tm_initialized = 0;

static TMStatus tm_grow_array(void **buffer, size_t *capacity, size_t elem_size, size_t min_capacity)
{
    void *new_buffer;
    size_t new_capacity;

    if (*capacity >= min_capacity) {
        return TM_OK;
    }

    new_capacity = (*capacity == 0) ? 8 : *capacity;
    while (new_capacity < min_capacity) {
        if (new_capacity > SIZE_MAX / 2) {
            return TM_ERR_ALLOC;
        }
        new_capacity *= 2;
    }

    if (elem_size != 0 && new_capacity > SIZE_MAX / elem_size) {
        return TM_ERR_ALLOC;
    }

    new_buffer = realloc(*buffer, new_capacity * elem_size);
    if (new_buffer == NULL) {
        return TM_ERR_ALLOC;
    }

    *buffer = new_buffer;
    *capacity = new_capacity;
    return TM_OK;
}

static void tm_swap_int(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

static void tm_swap_flag(unsigned char *a, unsigned char *b)
{
    unsigned char tmp = *a;
    *a = *b;
    *b = tmp;
}

static double tm_orient_value(const TMPoint *a, const TMPoint *b, const TMPoint *c)
{
    return orient2d((REAL *) a->xy, (REAL *) b->xy, (REAL *) c->xy);
}

static double tm_orient_coords(const TMPoint *a, const TMPoint *b, const double point[2])
{
    return orient2d((REAL *) a->xy, (REAL *) b->xy, (REAL *) point);
}

static double tm_incircle_value(const TMPoint *a, const TMPoint *b, const TMPoint *c, const TMPoint *d)
{
    return incircle((REAL *) a->xy, (REAL *) b->xy, (REAL *) c->xy, (REAL *) d->xy);
}

static int tm_point_nearly_on_segment(const TMPoint *a, const TMPoint *b, const double point[2])
{
    double dx = b->xy[0] - a->xy[0];
    double dy = b->xy[1] - a->xy[1];
    double len_sq = dx * dx + dy * dy;
    double scale;
    double tolerance;
    double t;
    double projection[2];
    double distance_sq;

    if (len_sq == 0.0) {
        return 0;
    }

    scale = fabs(a->xy[0]) + fabs(a->xy[1]) + fabs(b->xy[0]) + fabs(b->xy[1]) + fabs(point[0]) + fabs(point[1]);
    tolerance = 64.0 * DBL_EPSILON * ((scale > 1.0) ? scale : 1.0);

    t = ((point[0] - a->xy[0]) * dx + (point[1] - a->xy[1]) * dy) / len_sq;
    if (t < -tolerance || t > 1.0 + tolerance) {
        return 0;
    }

    projection[0] = a->xy[0] + t * dx;
    projection[1] = a->xy[1] + t * dy;
    distance_sq = (point[0] - projection[0]) * (point[0] - projection[0]) +
                  (point[1] - projection[1]) * (point[1] - projection[1]);
    return distance_sq <= tolerance * tolerance;
}

static int tm_point_nearly_equals_vertex(const TMPoint *vertex, const double point[2])
{
    double scale;
    double tolerance;
    double dx;
    double dy;
    double distance_sq;

    scale = fabs(vertex->xy[0]) + fabs(vertex->xy[1]) + fabs(point[0]) + fabs(point[1]);
    tolerance = 64.0 * DBL_EPSILON * ((scale > 1.0) ? scale : 1.0);
    dx = point[0] - vertex->xy[0];
    dy = point[1] - vertex->xy[1];
    distance_sq = dx * dx + dy * dy;
    return distance_sq <= tolerance * tolerance;
}

static int tm_compare_points_lex(const void *lhs, const void *rhs)
{
    const TMPoint *a = (const TMPoint *) lhs;
    const TMPoint *b = (const TMPoint *) rhs;

    if (a->xy[0] < b->xy[0]) {
        return -1;
    }
    if (a->xy[0] > b->xy[0]) {
        return 1;
    }
    if (a->xy[1] < b->xy[1]) {
        return -1;
    }
    if (a->xy[1] > b->xy[1]) {
        return 1;
    }
    if (a->original_index < b->original_index) {
        return -1;
    }
    if (a->original_index > b->original_index) {
        return 1;
    }
    return 0;
}

static int tm_points_equal(const TMPoint *a, const TMPoint *b)
{
    return a->xy[0] == b->xy[0] && a->xy[1] == b->xy[1];
}

static void tm_sort_three_ints(int values[3])
{
    if (values[0] > values[1]) {
        tm_swap_int(&values[0], &values[1]);
    }
    if (values[1] > values[2]) {
        tm_swap_int(&values[1], &values[2]);
    }
    if (values[0] > values[1]) {
        tm_swap_int(&values[0], &values[1]);
    }
}

static int tm_compare_triangles(const void *lhs, const void *rhs)
{
    const TMTriangle *a = (const TMTriangle *) lhs;
    const TMTriangle *b = (const TMTriangle *) rhs;
    int a_key[3] = {a->v[0], a->v[1], a->v[2]};
    int b_key[3] = {b->v[0], b->v[1], b->v[2]};
    int i;

    tm_sort_three_ints(a_key);
    tm_sort_three_ints(b_key);

    for (i = 0; i < 3; ++i) {
        if (a_key[i] < b_key[i]) {
            return -1;
        }
        if (a_key[i] > b_key[i]) {
            return 1;
        }
    }

    for (i = 0; i < 3; ++i) {
        if (a->v[i] < b->v[i]) {
            return -1;
        }
        if (a->v[i] > b->v[i]) {
            return 1;
        }
    }

    return 0;
}

static void tm_init_point(TMPoint *point, double x, double y, int original_index, TMVertexKind kind)
{
    point->xy[0] = x;
    point->xy[1] = y;
    point->original_index = original_index;
    point->kind = kind;
    point->incident_triangle = -1;
    point->protection_apex = -1;
    point->protection_side = TM_PROTECTION_SIDE_NONE;
    point->protection_level = 0;
}

static void tm_reset_triangle(TMTriangle *triangle)
{
    int edge;

    for (edge = 0; edge < 3; ++edge) {
        triangle->nbr[edge] = -1;
        triangle->constrained[edge] = 0;
    }
    triangle->dead = 0;
}

static void tm_normalize_triangle_ccw(TMTriangle *triangle, const TMPoint *points)
{
    double orientation = tm_orient_value(&points[triangle->v[0]], &points[triangle->v[1]], &points[triangle->v[2]]);

    if (orientation < 0.0) {
        tm_swap_int(&triangle->v[1], &triangle->v[2]);
        tm_swap_flag(&triangle->constrained[1], &triangle->constrained[2]);
    }
}

static int tm_make_triangle(TMTriangle *triangle, int a, int b, int c, const TMPoint *points)
{
    triangle->v[0] = a;
    triangle->v[1] = b;
    triangle->v[2] = c;
    tm_reset_triangle(triangle);
    tm_normalize_triangle_ccw(triangle, points);
    return tm_orient_value(&points[triangle->v[0]], &points[triangle->v[1]], &points[triangle->v[2]]) != 0.0;
}

static int tm_triangle_uses_supervertex(const TMWorkingTriangle *triangle, size_t real_point_count)
{
    return triangle->v[0] >= (int) real_point_count ||
           triangle->v[1] >= (int) real_point_count ||
           triangle->v[2] >= (int) real_point_count;
}

static void tm_compute_supertriangle(const TMPoint *points, size_t count, TMPoint supertriangle[3])
{
    double min_x = points[0].xy[0];
    double min_y = points[0].xy[1];
    double max_x = points[0].xy[0];
    double max_y = points[0].xy[1];
    double span;
    double mid_x;
    double mid_y;
    size_t i;

    for (i = 1; i < count; ++i) {
        if (points[i].xy[0] < min_x) {
            min_x = points[i].xy[0];
        }
        if (points[i].xy[0] > max_x) {
            max_x = points[i].xy[0];
        }
        if (points[i].xy[1] < min_y) {
            min_y = points[i].xy[1];
        }
        if (points[i].xy[1] > max_y) {
            max_y = points[i].xy[1];
        }
    }

    span = max_x - min_x;
    if (max_y - min_y > span) {
        span = max_y - min_y;
    }
    if (span == 0.0) {
        span = 1.0;
    }

    mid_x = 0.5 * (min_x + max_x);
    mid_y = 0.5 * (min_y + max_y);

    tm_init_point(&supertriangle[0], mid_x - 20.0 * span, mid_y - span, -1, TM_VERTEX_SUPER);
    tm_init_point(&supertriangle[1], mid_x, mid_y + 20.0 * span, -1, TM_VERTEX_SUPER);
    tm_init_point(&supertriangle[2], mid_x + 20.0 * span, mid_y - span, -1, TM_VERTEX_SUPER);

    if (tm_orient_value(&supertriangle[0], &supertriangle[1], &supertriangle[2]) < 0.0) {
        TMPoint tmp = supertriangle[1];
        supertriangle[1] = supertriangle[2];
        supertriangle[2] = tmp;
    }
}

static int tm_all_collinear(const TMPoint *points, size_t count)
{
    size_t i;

    for (i = 2; i < count; ++i) {
        if (tm_orient_value(&points[0], &points[1], &points[i]) != 0.0) {
            return 0;
        }
    }

    return 1;
}

static TMStatus tm_append_working_triangle(
    TMWorkingTriangle **triangles,
    size_t *triangle_count,
    size_t *triangle_capacity,
    int a,
    int b,
    int c,
    const TMPoint *points
)
{
    TMWorkingTriangle triangle;
    TMStatus status;

    triangle.v[0] = a;
    triangle.v[1] = b;
    triangle.v[2] = c;
    triangle.dead = 0;

    if (tm_orient_value(&points[a], &points[b], &points[c]) < 0.0) {
        tm_swap_int(&triangle.v[1], &triangle.v[2]);
    }

    if (tm_orient_value(&points[triangle.v[0]], &points[triangle.v[1]], &points[triangle.v[2]]) == 0.0) {
        return TM_OK;
    }

    status = tm_grow_array((void **) triangles, triangle_capacity, sizeof(**triangles), *triangle_count + 1);
    if (status != TM_OK) {
        return status;
    }

    (*triangles)[*triangle_count] = triangle;
    *triangle_count += 1;
    return TM_OK;
}

static TMStatus tm_add_boundary_edge(TMEdge **edges, size_t *edge_count, size_t *edge_capacity, int a, int b)
{
    TMStatus status;
    int lo;
    int hi;
    size_t i;

    lo = (a < b) ? a : b;
    hi = (a < b) ? b : a;

    for (i = 0; i < *edge_count; ++i) {
        if ((*edges)[i].active && (*edges)[i].a == lo && (*edges)[i].b == hi) {
            (*edges)[i].active = 0;
            return TM_OK;
        }
    }

    status = tm_grow_array((void **) edges, edge_capacity, sizeof(**edges), *edge_count + 1);
    if (status != TM_OK) {
        return status;
    }

    (*edges)[*edge_count].a = lo;
    (*edges)[*edge_count].b = hi;
    (*edges)[*edge_count].active = 1;
    *edge_count += 1;
    return TM_OK;
}

static TMStatus tm_append_point(TMMesh *mesh, const double point[2], int original_index, TMVertexKind kind, int *out_index)
{
    TMStatus status;

    status = tm_grow_array((void **) &mesh->points, &mesh->point_capacity, sizeof(*mesh->points), mesh->point_count + 1);
    if (status != TM_OK) {
        return status;
    }

    tm_init_point(&mesh->points[mesh->point_count], point[0], point[1], original_index, kind);
    if (out_index != NULL) {
        *out_index = (int) mesh->point_count;
    }
    mesh->point_count += 1;
    return TM_OK;
}

static TMStatus tm_append_triangle(TMMesh *mesh, const TMTriangle *triangle)
{
    TMStatus status;

    status = tm_grow_array(
        (void **) &mesh->triangles,
        &mesh->triangle_capacity,
        sizeof(*mesh->triangles),
        mesh->triangle_count + 1
    );
    if (status != TM_OK) {
        return status;
    }

    mesh->triangles[mesh->triangle_count] = *triangle;
    mesh->triangle_count += 1;
    return TM_OK;
}

static int tm_find_vertex_in_triangle(const TMTriangle *triangle, int vertex_index)
{
    int i;

    for (i = 0; i < 3; ++i) {
        if (triangle->v[i] == vertex_index) {
            return i;
        }
    }

    return -1;
}

static TMStatus tm_append_flip_candidate(
    TMFlipCandidate **candidates,
    size_t *count,
    size_t *capacity,
    int triangle_index,
    int edge
)
{
    TMStatus status;

    if (triangle_index < 0 || edge < 0 || edge > 2) {
        return TM_OK;
    }

    status = tm_grow_array((void **) candidates, capacity, sizeof(**candidates), *count + 1);
    if (status != TM_OK) {
        return status;
    }

    (*candidates)[*count].triangle_index = triangle_index;
    (*candidates)[*count].edge = edge;
    *count += 1;
    return TM_OK;
}

static TMStatus tm_legalize_point_star(TMMesh *mesh, int point_index, const int *seed_triangles, size_t seed_count)
{
    TMFlipCandidate *candidates = NULL;
    size_t candidate_count = 0;
    size_t candidate_capacity = 0;
    size_t i;
    TMStatus status = TM_OK;

    for (i = 0; i < seed_count; ++i) {
        int triangle_index = seed_triangles[i];
        int edge;

        if (triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count) {
            continue;
        }

        edge = tm_find_vertex_in_triangle(&mesh->triangles[triangle_index], point_index);
        if (edge < 0) {
            continue;
        }

        status = tm_append_flip_candidate(
            &candidates,
            &candidate_count,
            &candidate_capacity,
            triangle_index,
            edge
        );
        if (status != TM_OK) {
            free(candidates);
            return status;
        }
    }

    while (candidate_count != 0) {
        TMFlipCandidate candidate = candidates[candidate_count - 1];
        int triangle_index = candidate.triangle_index;
        int edge = candidate.edge;
        int neighbor_index;

        candidate_count -= 1;

        if (triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count) {
            continue;
        }

        edge = tm_find_vertex_in_triangle(&mesh->triangles[triangle_index], point_index);
        if (edge < 0) {
            continue;
        }

        neighbor_index = mesh->triangles[triangle_index].nbr[edge];
        if (neighbor_index < 0) {
            continue;
        }
        if (mesh->triangles[triangle_index].constrained[edge]) {
            continue;
        }
        if (tm_edge_is_locally_delaunay(mesh, triangle_index, edge)) {
            continue;
        }

        status = tm_flip_edge(mesh, triangle_index, edge);
        if (status != TM_OK) {
            free(candidates);
            return status;
        }

        status = tm_append_flip_candidate(
            &candidates,
            &candidate_count,
            &candidate_capacity,
            triangle_index,
            tm_find_vertex_in_triangle(&mesh->triangles[triangle_index], point_index)
        );
        if (status == TM_OK) {
            status = tm_append_flip_candidate(
                &candidates,
                &candidate_count,
                &candidate_capacity,
                neighbor_index,
                tm_find_vertex_in_triangle(&mesh->triangles[neighbor_index], point_index)
            );
        }
        if (status != TM_OK) {
            free(candidates);
            return status;
        }
    }

    free(candidates);
    return TM_OK;
}

static TMStatus tm_set_triangle_edge_neighbor(TMTriangle *triangle, int a, int b, int neighbor)
{
    int edge = tm_find_edge_in_triangle(triangle, a, b);

    if (edge < 0) {
        return TM_ERR_INVALID_MESH;
    }

    triangle->nbr[edge] = neighbor;
    return TM_OK;
}

static TMStatus tm_set_triangle_edge_constraint(TMTriangle *triangle, int a, int b, unsigned char constrained)
{
    int edge = tm_find_edge_in_triangle(triangle, a, b);

    if (edge < 0) {
        return TM_ERR_INVALID_MESH;
    }

    triangle->constrained[edge] = constrained;
    return TM_OK;
}

static TMStatus tm_copy_triangle_edge_state(const TMTriangle *source, int a, int b, TMTriangle *target)
{
    int source_edge = tm_find_edge_in_triangle(source, a, b);
    int target_edge = tm_find_edge_in_triangle(target, a, b);

    if (source_edge < 0 || target_edge < 0) {
        return TM_ERR_INVALID_MESH;
    }

    target->nbr[target_edge] = source->nbr[source_edge];
    target->constrained[target_edge] = source->constrained[source_edge];
    return TM_OK;
}

static TMStatus tm_redirect_neighbor_reference(
    TMMesh *mesh,
    int triangle_index,
    int a,
    int b,
    int old_neighbor,
    int new_neighbor
)
{
    int edge;

    if (triangle_index < 0) {
        return TM_OK;
    }
    if ((size_t) triangle_index >= mesh->triangle_count) {
        return TM_ERR_INVALID_MESH;
    }
    if (old_neighbor == new_neighbor) {
        return TM_OK;
    }

    edge = tm_find_edge_in_triangle(&mesh->triangles[triangle_index], a, b);
    if (edge < 0 || mesh->triangles[triangle_index].nbr[edge] != old_neighbor) {
        return TM_ERR_INVALID_MESH;
    }

    mesh->triangles[triangle_index].nbr[edge] = new_neighbor;
    return TM_OK;
}

static int tm_edge_matches_live_segment(const TMMesh *mesh, int a, int b)
{
    size_t segment_index;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;

    if (mesh == NULL) {
        return 0;
    }

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &mesh->segments[segment_index];
        int seg_lo;
        int seg_hi;

        if (!segment->live) {
            continue;
        }

        seg_lo = (segment->v[0] < segment->v[1]) ? segment->v[0] : segment->v[1];
        seg_hi = (segment->v[0] < segment->v[1]) ? segment->v[1] : segment->v[0];
        if (seg_lo == lo && seg_hi == hi) {
            return 1;
        }
    }

    return 0;
}

static TMStatus tm_link_triangles_across_edge(TMMesh *mesh, int first_index, int second_index, int a, int b)
{
    TMStatus status;

    if (first_index < 0 || second_index < 0 ||
        (size_t) first_index >= mesh->triangle_count ||
        (size_t) second_index >= mesh->triangle_count) {
        return TM_ERR_INVALID_MESH;
    }

    status = tm_set_triangle_edge_neighbor(&mesh->triangles[first_index], a, b, second_index);
    if (status != TM_OK) {
        return status;
    }

    return tm_set_triangle_edge_neighbor(&mesh->triangles[second_index], a, b, first_index);
}

static void tm_assign_incident_triangle(TMMesh *mesh, int point_index, int triangle_index)
{
    if (point_index < 0 || triangle_index < 0) {
        return;
    }
    if ((size_t) point_index >= mesh->point_count || (size_t) triangle_index >= mesh->triangle_count) {
        return;
    }

    mesh->points[point_index].incident_triangle = triangle_index;
}

static void tm_compact_triangles(TMMesh *mesh)
{
    size_t read_index;
    size_t write_index = 0;

    for (read_index = 0; read_index < mesh->triangle_count; ++read_index) {
        if (!mesh->triangles[read_index].dead) {
            if (write_index != read_index) {
                mesh->triangles[write_index] = mesh->triangles[read_index];
            }
            mesh->triangles[write_index].dead = 0;
            write_index += 1;
        }
    }

    mesh->triangle_count = write_index;
}

static TMStatus tm_append_edge_pair(TMEdgePair **pairs, size_t *pair_count, size_t *pair_capacity, int a, int b)
{
    TMStatus status;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;
    size_t i;

    for (i = 0; i < *pair_count; ++i) {
        if ((*pairs)[i].a == lo && (*pairs)[i].b == hi) {
            return TM_OK;
        }
    }

    status = tm_grow_array((void **) pairs, pair_capacity, sizeof(**pairs), *pair_count + 1);
    if (status != TM_OK) {
        return status;
    }

    (*pairs)[*pair_count].a = lo;
    (*pairs)[*pair_count].b = hi;
    *pair_count += 1;
    return TM_OK;
}

static void tm_remove_edge_pair(TMEdgePair *pairs, size_t *pair_count, int a, int b)
{
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;
    size_t i;

    for (i = 0; i < *pair_count; ++i) {
        if (pairs[i].a == lo && pairs[i].b == hi) {
            pairs[i] = pairs[*pair_count - 1];
            *pair_count -= 1;
            return;
        }
    }
}

static TMStatus tm_collect_constrained_edge_pairs(const TMMesh *mesh, TMEdgePair **out_pairs, size_t *out_pair_count)
{
    TMEdgePair *pairs = NULL;
    size_t pair_count = 0;
    size_t pair_capacity = 0;
    size_t tri_index;
    TMStatus status = TM_OK;

    *out_pairs = NULL;
    *out_pair_count = 0;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            int a;
            int b;

            if (!mesh->triangles[tri_index].constrained[edge]) {
                continue;
            }

            tm_triangle_edge_vertices(&mesh->triangles[tri_index], edge, &a, &b);
            status = tm_append_edge_pair(&pairs, &pair_count, &pair_capacity, a, b);
            if (status != TM_OK) {
                free(pairs);
                return status;
            }
        }
    }

    *out_pairs = pairs;
    *out_pair_count = pair_count;
    return TM_OK;
}

static void tm_clear_constraints(TMMesh *mesh)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            mesh->triangles[tri_index].constrained[edge] = 0;
        }
    }
}

static TMStatus tm_apply_constraint_pairs(TMMesh *mesh, const TMEdgePair *pairs, size_t pair_count)
{
    size_t pair_index;

    tm_clear_constraints(mesh);

    for (pair_index = 0; pair_index < pair_count; ++pair_index) {
        size_t tri_index;
        int found = 0;

        for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
            int edge = tm_find_edge_in_triangle(&mesh->triangles[tri_index], pairs[pair_index].a, pairs[pair_index].b);

            if (edge >= 0) {
                mesh->triangles[tri_index].constrained[edge] = 1;
                found = 1;
            }
        }

        if (!found) {
            return TM_ERR_INVALID_MESH;
        }
    }

    return TM_OK;
}

static int tm_point_inside_triangle(const TMMesh *mesh, const TMTriangle *triangle, const double point[2], int *out_edge)
{
    int edge;
    int on_edge = -1;

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double side;

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        side = tm_orient_coords(&mesh->points[a], &mesh->points[b], point);
        if (side < 0.0) {
            return 0;
        }
        if (side == 0.0) {
            on_edge = edge;
        }
    }

    if (out_edge != NULL) {
        *out_edge = on_edge;
    }
    return 1;
}

void tm_initialize(void)
{
    if (!tm_initialized) {
        exactinit();
        tm_initialized = 1;
    }
}

const char *tm_internal_status_string(TMStatus status)
{
    switch (status) {
        case TM_OK:
            return "ok";
        case TM_ERR_ALLOC:
            return "allocation failed";
        case TM_ERR_IO:
            return "I/O error";
        case TM_ERR_PARSE:
            return "input parse error";
        case TM_ERR_DUPLICATE:
            return "duplicate points are not allowed";
        case TM_ERR_TOO_FEW_POINTS:
            return "at least three unique points are required";
        case TM_ERR_COLLINEAR:
            return "all points are collinear";
        case TM_ERR_INVALID_PSLG:
            return "invalid PSLG input";
        case TM_ERR_INVALID_MESH:
            return "invalid mesh topology";
        case TM_ERR_INTERNAL:
            return "internal triangulation error";
    }

    return "unknown error";
}

TMStatus tm_read_points_file(const char *path, TMPoint **out_points, size_t *out_count)
{
    FILE *stream;
    TMPoint *points = NULL;
    size_t point_count = 0;
    size_t point_capacity = 0;
    char line[4096];
    TMStatus status = TM_OK;

    if (out_points == NULL || out_count == NULL) {
        return TM_ERR_INTERNAL;
    }

    *out_points = NULL;
    *out_count = 0;

    stream = fopen(path, "r");
    if (stream == NULL) {
        return TM_ERR_IO;
    }

    while (fgets(line, sizeof(line), stream) != NULL) {
        char *cursor = line;
        double x;
        double y;
        char extra;
        int matched;

        while (*cursor != '\0' && isspace((unsigned char) *cursor)) {
            cursor++;
        }

        if (*cursor == '\0' || *cursor == '\n' || *cursor == '#') {
            continue;
        }

        matched = sscanf(cursor, "%lf %lf %c", &x, &y, &extra);
        if (matched != 2 || !isfinite(x) || !isfinite(y)) {
            status = TM_ERR_PARSE;
            break;
        }

        status = tm_grow_array((void **) &points, &point_capacity, sizeof(*points), point_count + 1);
        if (status != TM_OK) {
            break;
        }

        tm_init_point(&points[point_count], x, y, (int) point_count, TM_VERTEX_INPUT);
        point_count += 1;
    }

    if (status == TM_OK && ferror(stream)) {
        status = TM_ERR_IO;
    }

    fclose(stream);

    if (status != TM_OK) {
        free(points);
        return status;
    }

    *out_points = points;
    *out_count = point_count;
    return TM_OK;
}

TMStatus tm_build_mesh(const TMPoint *input_points, size_t input_count, TMMesh *out_mesh)
{
    TMPoint *working_points = NULL;
    TMWorkingTriangle *triangles = NULL;
    size_t triangle_count = 0;
    size_t triangle_capacity = 0;
    size_t final_count = 0;
    size_t i;
    TMStatus status = TM_OK;

    if (out_mesh == NULL || input_points == NULL) {
        return TM_ERR_INTERNAL;
    }

    memset(out_mesh, 0, sizeof(*out_mesh));

    if (input_count < 3) {
        return TM_ERR_TOO_FEW_POINTS;
    }

    tm_initialize();

    working_points = (TMPoint *) malloc((input_count + 3) * sizeof(*working_points));
    if (working_points == NULL) {
        return TM_ERR_ALLOC;
    }

    memcpy(working_points, input_points, input_count * sizeof(*working_points));
    qsort(working_points, input_count, sizeof(*working_points), tm_compare_points_lex);

    for (i = 1; i < input_count; ++i) {
        if (tm_points_equal(&working_points[i - 1], &working_points[i])) {
            status = TM_ERR_DUPLICATE;
            goto cleanup;
        }
    }

    if (tm_all_collinear(working_points, input_count)) {
        status = TM_ERR_COLLINEAR;
        goto cleanup;
    }

    tm_compute_supertriangle(working_points, input_count, &working_points[input_count]);

    status = tm_append_working_triangle(
        &triangles,
        &triangle_count,
        &triangle_capacity,
        (int) input_count,
        (int) input_count + 1,
        (int) input_count + 2,
        working_points
    );
    if (status != TM_OK) {
        goto cleanup;
    }

    for (i = 0; i < input_count; ++i) {
        TMEdge *edges = NULL;
        size_t edge_count = 0;
        size_t edge_capacity = 0;
        size_t tri_index;
        size_t bad_count = 0;

        for (tri_index = 0; tri_index < triangle_count; ++tri_index) {
            TMWorkingTriangle *triangle = &triangles[tri_index];
            double in_value;

            if (triangle->dead) {
                continue;
            }

            in_value = tm_incircle_value(
                &working_points[triangle->v[0]],
                &working_points[triangle->v[1]],
                &working_points[triangle->v[2]],
                &working_points[i]
            );

            if (in_value >= 0.0) {
                triangle->dead = 1;
                bad_count += 1;

                status = tm_add_boundary_edge(&edges, &edge_count, &edge_capacity, triangle->v[0], triangle->v[1]);
                if (status != TM_OK) {
                    free(edges);
                    goto cleanup;
                }
                status = tm_add_boundary_edge(&edges, &edge_count, &edge_capacity, triangle->v[1], triangle->v[2]);
                if (status != TM_OK) {
                    free(edges);
                    goto cleanup;
                }
                status = tm_add_boundary_edge(&edges, &edge_count, &edge_capacity, triangle->v[2], triangle->v[0]);
                if (status != TM_OK) {
                    free(edges);
                    goto cleanup;
                }
            }
        }

        if (bad_count == 0) {
            free(edges);
            status = TM_ERR_INTERNAL;
            goto cleanup;
        }

        for (tri_index = 0; tri_index < edge_count; ++tri_index) {
            if (!edges[tri_index].active) {
                continue;
            }

            status = tm_append_working_triangle(
                &triangles,
                &triangle_count,
                &triangle_capacity,
                edges[tri_index].a,
                edges[tri_index].b,
                (int) i,
                working_points
            );
            if (status != TM_OK) {
                free(edges);
                goto cleanup;
            }
        }

        free(edges);
    }

    for (i = 0; i < triangle_count; ++i) {
        if (!triangles[i].dead && !tm_triangle_uses_supervertex(&triangles[i], input_count)) {
            final_count += 1;
        }
    }

    if (final_count == 0) {
        status = TM_ERR_INTERNAL;
        goto cleanup;
    }

    out_mesh->points = (TMPoint *) malloc(input_count * sizeof(*out_mesh->points));
    if (out_mesh->points == NULL) {
        status = TM_ERR_ALLOC;
        goto cleanup;
    }
    out_mesh->point_capacity = input_count;
    out_mesh->point_count = input_count;
    out_mesh->segments = NULL;
    out_mesh->segment_count = 0;
    out_mesh->segment_capacity = 0;

    for (i = 0; i < input_count; ++i) {
        out_mesh->points[i] = input_points[i];
        out_mesh->points[i].kind = TM_VERTEX_INPUT;
        out_mesh->points[i].incident_triangle = -1;
        out_mesh->points[i].protection_apex = -1;
        out_mesh->points[i].protection_side = TM_PROTECTION_SIDE_NONE;
        out_mesh->points[i].protection_level = 0;
    }

    out_mesh->triangles = (TMTriangle *) malloc(final_count * sizeof(*out_mesh->triangles));
    if (out_mesh->triangles == NULL) {
        status = TM_ERR_ALLOC;
        goto cleanup;
    }
    out_mesh->triangle_capacity = final_count;
    out_mesh->triangle_count = 0;

    for (i = 0; i < triangle_count; ++i) {
        if (!triangles[i].dead && !tm_triangle_uses_supervertex(&triangles[i], input_count)) {
            TMTriangle triangle;

            if (!tm_make_triangle(
                    &triangle,
                    working_points[triangles[i].v[0]].original_index,
                    working_points[triangles[i].v[1]].original_index,
                    working_points[triangles[i].v[2]].original_index,
                    out_mesh->points
                )) {
                status = TM_ERR_INTERNAL;
                goto cleanup;
            }

            out_mesh->triangles[out_mesh->triangle_count] = triangle;
            out_mesh->triangle_count += 1;
        }
    }

    qsort(out_mesh->triangles, out_mesh->triangle_count, sizeof(*out_mesh->triangles), tm_compare_triangles);
    status = tm_rebuild_topology(out_mesh);
    if (status != TM_OK) {
        goto cleanup;
    }

cleanup:
    free(triangles);
    free(working_points);

    if (status != TM_OK) {
        tm_free_mesh(out_mesh);
    }

    return status;
}

TMStatus tm_rebuild_topology(TMMesh *mesh)
{
    size_t tri_index;

    if (mesh == NULL || mesh->points == NULL || mesh->triangles == NULL) {
        return TM_ERR_INTERNAL;
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        tm_normalize_triangle_ccw(&mesh->triangles[tri_index], mesh->points);
        mesh->triangles[tri_index].dead = 0;
        for (edge = 0; edge < 3; ++edge) {
            mesh->triangles[tri_index].nbr[edge] = -1;
        }
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        size_t other_index;

        for (other_index = tri_index + 1; other_index < mesh->triangle_count; ++other_index) {
            int edge;

            for (edge = 0; edge < 3; ++edge) {
                int a;
                int b;
                int other_edge;

                if (mesh->triangles[tri_index].nbr[edge] != -1) {
                    continue;
                }

                tm_triangle_edge_vertices(&mesh->triangles[tri_index], edge, &a, &b);
                other_edge = tm_find_edge_in_triangle(&mesh->triangles[other_index], a, b);
                if (other_edge >= 0) {
                    mesh->triangles[tri_index].nbr[edge] = (int) other_index;
                    mesh->triangles[other_index].nbr[other_edge] = (int) tri_index;
                }
            }
        }
    }

    for (tri_index = 0; tri_index < mesh->point_count; ++tri_index) {
        mesh->points[tri_index].incident_triangle = -1;
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int vertex;

        for (vertex = 0; vertex < 3; ++vertex) {
            int point_index = mesh->triangles[tri_index].v[vertex];

            if (point_index < 0 || (size_t) point_index >= mesh->point_count) {
                return TM_ERR_INVALID_MESH;
            }
            if (mesh->points[point_index].incident_triangle == -1) {
                mesh->points[point_index].incident_triangle = (int) tri_index;
            }
        }
    }

    return TM_OK;
}

void tm_triangle_edge_vertices(const TMTriangle *triangle, int edge, int *out_a, int *out_b)
{
    if (out_a != NULL) {
        *out_a = triangle->v[(edge + 1) % 3];
    }
    if (out_b != NULL) {
        *out_b = triangle->v[(edge + 2) % 3];
    }
}

int tm_find_edge_in_triangle(const TMTriangle *triangle, int a, int b)
{
    int edge;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;

    for (edge = 0; edge < 3; ++edge) {
        int edge_a;
        int edge_b;
        int edge_lo;
        int edge_hi;

        tm_triangle_edge_vertices(triangle, edge, &edge_a, &edge_b);
        edge_lo = (edge_a < edge_b) ? edge_a : edge_b;
        edge_hi = (edge_a < edge_b) ? edge_b : edge_a;

        if (edge_lo == lo && edge_hi == hi) {
            return edge;
        }
    }

    return -1;
}

int tm_triangle_contains_vertex(const TMTriangle *triangle, int vertex_index)
{
    return triangle->v[0] == vertex_index || triangle->v[1] == vertex_index || triangle->v[2] == vertex_index;
}

int tm_edge_is_flippable(const TMMesh *mesh, int triangle_index, int edge)
{
    const TMTriangle *triangle;
    const TMTriangle *neighbor;
    int a;
    int b;
    int neighbor_index;
    int neighbor_edge;
    int c;
    int d;

    if (mesh == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count || edge < 0 || edge > 2) {
        return 0;
    }

    triangle = &mesh->triangles[triangle_index];
    neighbor_index = triangle->nbr[edge];
    if (neighbor_index < 0 || triangle->constrained[edge]) {
        return 0;
    }

    neighbor = &mesh->triangles[neighbor_index];
    tm_triangle_edge_vertices(triangle, edge, &a, &b);
    if (tm_edge_matches_live_segment(mesh, a, b)) {
        return 0;
    }
    neighbor_edge = tm_find_edge_in_triangle(neighbor, a, b);
    if (neighbor_edge < 0 || neighbor->constrained[neighbor_edge]) {
        return 0;
    }

    c = triangle->v[edge];
    d = neighbor->v[neighbor_edge];
    return tm_orient_value(&mesh->points[c], &mesh->points[a], &mesh->points[d]) > 0.0 &&
           tm_orient_value(&mesh->points[c], &mesh->points[d], &mesh->points[b]) > 0.0;
}

int tm_edge_is_locally_delaunay(const TMMesh *mesh, int triangle_index, int edge)
{
    const TMTriangle *triangle;
    const TMTriangle *neighbor;
    int a;
    int b;
    int neighbor_index;
    int neighbor_edge;
    int c;
    int d;

    if (mesh == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count || edge < 0 || edge > 2) {
        return 0;
    }

    if (!tm_edge_is_flippable(mesh, triangle_index, edge)) {
        return 1;
    }

    triangle = &mesh->triangles[triangle_index];
    neighbor_index = triangle->nbr[edge];
    neighbor = &mesh->triangles[neighbor_index];
    tm_triangle_edge_vertices(triangle, edge, &a, &b);
    neighbor_edge = tm_find_edge_in_triangle(neighbor, a, b);
    c = triangle->v[edge];
    d = neighbor->v[neighbor_edge];
    return tm_incircle_value(&mesh->points[a], &mesh->points[b], &mesh->points[c], &mesh->points[d]) <= 0.0;
}

static int tm_find_live_edge_segment(const TMMesh *mesh, int a, int b, size_t *out_segment_index)
{
    size_t segment_index;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        int sa;
        int sb;
        int slo;
        int shi;

        if (!mesh->segments[segment_index].live) {
            continue;
        }

        sa = mesh->segments[segment_index].v[0];
        sb = mesh->segments[segment_index].v[1];
        slo = (sa < sb) ? sa : sb;
        shi = (sa < sb) ? sb : sa;
        if (slo == lo && shi == hi) {
            if (out_segment_index != NULL) {
                *out_segment_index = segment_index;
            }
            return 1;
        }
    }

    return 0;
}

static int tm_shell_edge_is_exempt(const TMMesh *mesh, int a, int b)
{
    const TMPoint *pa = &mesh->points[a];
    const TMPoint *pb = &mesh->points[b];

    if (pa->protection_apex < 0 || pb->protection_apex < 0) {
        return 0;
    }
    if (pa->protection_apex != pb->protection_apex) {
        return 0;
    }
    if (pa->protection_side == TM_PROTECTION_SIDE_NONE ||
        pb->protection_side == TM_PROTECTION_SIDE_NONE ||
        pa->protection_side == pb->protection_side) {
        return 0;
    }
    if (pa->protection_level == 0 || pb->protection_level == 0) {
        return 0;
    }
    return 1;
}

static int tm_shell_chain_edge_is_exempt(const TMMesh *mesh, int a, int b)
{
    const TMPoint *pa = &mesh->points[a];
    const TMPoint *pb = &mesh->points[b];
    unsigned int level_lo;
    unsigned int level_hi;
    size_t segment_index;

    if (pa->protection_apex < 0 || pb->protection_apex < 0) {
        return 0;
    }
    if (pa->protection_apex != pb->protection_apex) {
        return 0;
    }
    if (pa->protection_side == TM_PROTECTION_SIDE_NONE ||
        pb->protection_side == TM_PROTECTION_SIDE_NONE ||
        pa->protection_side != pb->protection_side) {
        return 0;
    }
    if (pa->protection_level == 0 || pb->protection_level == 0) {
        return 0;
    }

    level_lo = (pa->protection_level < pb->protection_level) ? pa->protection_level : pb->protection_level;
    level_hi = (pa->protection_level < pb->protection_level) ? pb->protection_level : pa->protection_level;
    if (level_hi - level_lo != 1) {
        return 0;
    }
    if (!tm_find_live_edge_segment(mesh, a, b, &segment_index)) {
        return 0;
    }

    return mesh->segments[segment_index].is_protected != 0;
}

static int tm_triangle_shortest_edge_is_exempt(const TMMesh *mesh, const TMTriangle *triangle)
{
    int edge;
    int best_a = -1;
    int best_b = -1;
    double best_length_sq = 0.0;

    if (mesh == NULL || triangle == NULL) {
        return 0;
    }

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double dx;
        double dy;
        double length_sq;

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        dx = mesh->points[b].xy[0] - mesh->points[a].xy[0];
        dy = mesh->points[b].xy[1] - mesh->points[a].xy[1];
        length_sq = dx * dx + dy * dy;

        if (edge == 0 || length_sq < best_length_sq) {
            best_length_sq = length_sq;
            best_a = a;
            best_b = b;
        }
    }

    return tm_shell_edge_is_exempt(mesh, best_a, best_b);
}

int tm_triangle_is_quality_exempt(const TMMesh *mesh, int triangle_index)
{
    const TMTriangle *triangle;
    int edge;

    if (mesh == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count) {
        return 0;
    }

    triangle = &mesh->triangles[triangle_index];
    if (tm_triangle_shortest_edge_is_exempt(mesh, triangle)) {
        return 1;
    }

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        if (tm_shell_chain_edge_is_exempt(mesh, a, b)) {
            return 1;
        }
    }

    return 0;
}

size_t tm_count_protected_corners(const TMMesh *mesh)
{
    size_t point_index;
    size_t protected_corner_count = 0;

    if (mesh == NULL) {
        return 0;
    }

    for (point_index = 0; point_index < mesh->point_count; ++point_index) {
        size_t segment_index;
        size_t count = 0;

        for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
            if (!mesh->segments[segment_index].live || !mesh->segments[segment_index].is_protected) {
                continue;
            }
            if (mesh->segments[segment_index].protected_apex == (int) point_index) {
                count += 1;
            }
        }

        if (count >= 2) {
            protected_corner_count += 1;
        }
    }

    return protected_corner_count;
}

TMStatus tm_locate_point(const TMMesh *mesh, const double point[2], int start_triangle, TMLocation *out_location)
{
    int current;
    int steps;
    int max_steps;

    if (mesh == NULL || point == NULL || out_location == NULL || mesh->triangle_count == 0) {
        return TM_ERR_INTERNAL;
    }

    current = (start_triangle >= 0 && (size_t) start_triangle < mesh->triangle_count) ? start_triangle : 0;
    max_steps = (int) (mesh->triangle_count * 4 + 4);

    for (steps = 0; steps < max_steps; ++steps) {
        const TMTriangle *triangle = &mesh->triangles[current];
        int edge;
        int boundary_edge = -1;
        int on_edge = -1;

        for (edge = 0; edge < 3; ++edge) {
            int a;
            int b;
            double side;
            int near_edge;

            tm_triangle_edge_vertices(triangle, edge, &a, &b);
            side = tm_orient_coords(&mesh->points[a], &mesh->points[b], point);
            near_edge = tm_point_nearly_on_segment(&mesh->points[a], &mesh->points[b], point);
            if (side < 0.0 && !near_edge) {
                boundary_edge = edge;
                break;
            }
            if (side == 0.0 || near_edge) {
                on_edge = edge;
            }
        }

        if (boundary_edge == -1) {
            out_location->triangle = current;
            out_location->edge = on_edge;
            out_location->on_edge = (on_edge >= 0);
            return TM_OK;
        }

        if (triangle->nbr[boundary_edge] < 0) {
            out_location->triangle = current;
            out_location->edge = boundary_edge;
            out_location->on_edge = 0;
            return TM_OK;
        }

        current = triangle->nbr[boundary_edge];
    }

    return TM_ERR_INVALID_MESH;
}

TMStatus tm_flip_edge(TMMesh *mesh, int triangle_index, int edge)
{
    TMTriangle first;
    TMTriangle second;
    TMTriangle triangle;
    TMTriangle neighbor;
    int neighbor_index;
    int neighbor_edge;
    int a;
    int b;
    int c;
    int d;
    int neighbor_ca;
    int neighbor_bc;
    int neighbor_db;
    int neighbor_ad;
    unsigned char constrained_ca;
    unsigned char constrained_bc;
    unsigned char constrained_db;
    unsigned char constrained_ad;
    TMStatus status;

    if (mesh == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count || edge < 0 || edge > 2) {
        tm_set_pslg_error_detail("flip edge: invalid triangle=%d edge=%d", triangle_index, edge);
        return TM_ERR_INTERNAL;
    }

    neighbor_index = mesh->triangles[triangle_index].nbr[edge];
    if (neighbor_index < 0 || mesh->triangles[triangle_index].constrained[edge]) {
        tm_set_pslg_error_detail(
            "flip edge: triangle %d edge %d has no flippable neighbor (neighbor=%d constrained=%d)",
            triangle_index,
            edge,
            neighbor_index,
            (int) mesh->triangles[triangle_index].constrained[edge]
        );
        return TM_ERR_INTERNAL;
    }

    tm_triangle_edge_vertices(&mesh->triangles[triangle_index], edge, &a, &b);
    if (tm_edge_matches_live_segment(mesh, a, b)) {
        tm_set_pslg_error_detail(
            "flip edge: edge %d-%d is a live constrained segment",
            a,
            b
        );
        return TM_ERR_INTERNAL;
    }
    neighbor_edge = tm_find_edge_in_triangle(&mesh->triangles[neighbor_index], a, b);
    if (neighbor_edge < 0 || mesh->triangles[neighbor_index].constrained[neighbor_edge]) {
        tm_set_pslg_error_detail(
            "flip edge: neighbor %d does not mirror edge %d-%d (neighbor_edge=%d constrained=%d)",
            neighbor_index,
            a,
            b,
            neighbor_edge,
            neighbor_edge >= 0 ? (int) mesh->triangles[neighbor_index].constrained[neighbor_edge] : -1
        );
        return TM_ERR_INTERNAL;
    }

    c = mesh->triangles[triangle_index].v[edge];
    d = mesh->triangles[neighbor_index].v[neighbor_edge];
    triangle = mesh->triangles[triangle_index];
    neighbor = mesh->triangles[neighbor_index];
    neighbor_ca = triangle.nbr[tm_find_edge_in_triangle(&triangle, c, a)];
    neighbor_bc = triangle.nbr[tm_find_edge_in_triangle(&triangle, b, c)];
    neighbor_db = neighbor.nbr[tm_find_edge_in_triangle(&neighbor, d, b)];
    neighbor_ad = neighbor.nbr[tm_find_edge_in_triangle(&neighbor, a, d)];
    constrained_ca = triangle.constrained[tm_find_edge_in_triangle(&triangle, c, a)];
    constrained_bc = triangle.constrained[tm_find_edge_in_triangle(&triangle, b, c)];
    constrained_db = neighbor.constrained[tm_find_edge_in_triangle(&neighbor, d, b)];
    constrained_ad = neighbor.constrained[tm_find_edge_in_triangle(&neighbor, a, d)];

    if (tm_orient_value(&mesh->points[c], &mesh->points[a], &mesh->points[d]) <= 0.0 ||
        tm_orient_value(&mesh->points[c], &mesh->points[d], &mesh->points[b]) <= 0.0) {
        tm_set_pslg_error_detail(
            "flip edge: non-convex quad for diagonal %d-%d with opposite vertices %d and %d",
            a,
            b,
            c,
            d
        );
        return TM_ERR_INTERNAL;
    }

    if (!tm_make_triangle(&first, c, a, d, mesh->points) || !tm_make_triangle(&second, c, d, b, mesh->points)) {
        tm_set_pslg_error_detail(
            "flip edge: replacement triangles degenerate for diagonal %d-%d with opposite vertices %d and %d",
            a,
            b,
            c,
            d
        );
        return TM_ERR_INTERNAL;
    }

    status = tm_copy_triangle_edge_state(&triangle, c, a, &first);
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&neighbor, a, d, &first);
    }
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&triangle, b, c, &second);
    }
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&neighbor, d, b, &second);
    }
    if (status == TM_OK) {
        mesh->triangles[triangle_index] = first;
        mesh->triangles[neighbor_index] = second;
        status = tm_link_triangles_across_edge(mesh, triangle_index, neighbor_index, c, d);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[triangle_index], c, a, constrained_ca);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[triangle_index], a, d, constrained_ad);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[neighbor_index], b, c, constrained_bc);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[neighbor_index], d, b, constrained_db);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_neighbor(&mesh->triangles[triangle_index], c, a, neighbor_ca);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_neighbor(&mesh->triangles[triangle_index], a, d, neighbor_ad);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_neighbor(&mesh->triangles[neighbor_index], b, c, neighbor_bc);
    }
    if (status == TM_OK) {
        status = tm_set_triangle_edge_neighbor(&mesh->triangles[neighbor_index], d, b, neighbor_db);
    }
    if (status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, triangle_index, neighbor_index, c, d);
    }
    if (status == TM_OK) {
        status = tm_redirect_neighbor_reference(mesh, neighbor_bc, b, c, triangle_index, neighbor_index);
    }
    if (status == TM_OK) {
        status = tm_redirect_neighbor_reference(mesh, neighbor_ad, a, d, neighbor_index, triangle_index);
    }
    if (status != TM_OK) {
        tm_set_pslg_error_detail(
            "flip edge: local topology update failed after flipping diagonal %d-%d",
            a,
            b
        );
        return status;
    }

    tm_assign_incident_triangle(mesh, a, triangle_index);
    tm_assign_incident_triangle(mesh, b, neighbor_index);
    tm_assign_incident_triangle(mesh, c, triangle_index);
    tm_assign_incident_triangle(mesh, d, neighbor_index);
    return TM_OK;
}

TMStatus tm_insert_point_in_triangle(
    TMMesh *mesh,
    int triangle_index,
    const double point[2],
    TMVertexKind kind,
    int *out_point_index
)
{
    TMTriangle first;
    TMTriangle second;
    TMTriangle third;
    TMTriangle triangle;
    int v0;
    int v1;
    int v2;
    int point_index;
    int second_index;
    int third_index;
    int seed_triangles[3];
    TMStatus status;
    int on_edge = -1;

    if (mesh == NULL || point == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count) {
        return TM_ERR_INTERNAL;
    }

    if (!tm_point_inside_triangle(mesh, &mesh->triangles[triangle_index], point, &on_edge) || on_edge >= 0) {
        return TM_ERR_INTERNAL;
    }

    v0 = mesh->triangles[triangle_index].v[0];
    v1 = mesh->triangles[triangle_index].v[1];
    v2 = mesh->triangles[triangle_index].v[2];
    if (kind == TM_VERTEX_TRIANGLE_SPLIT) {
        if (tm_point_nearly_equals_vertex(&mesh->points[v0], point)) {
            if (out_point_index != NULL) {
                *out_point_index = v0;
            }
            return TM_OK;
        }
        if (tm_point_nearly_equals_vertex(&mesh->points[v1], point)) {
            if (out_point_index != NULL) {
                *out_point_index = v1;
            }
            return TM_OK;
        }
        if (tm_point_nearly_equals_vertex(&mesh->points[v2], point)) {
            if (out_point_index != NULL) {
                *out_point_index = v2;
            }
            return TM_OK;
        }
    }
    triangle = mesh->triangles[triangle_index];
    second_index = (int) mesh->triangle_count;
    third_index = second_index + 1;

    status = tm_grow_array((void **) &mesh->points, &mesh->point_capacity, sizeof(*mesh->points), mesh->point_count + 1);
    if (status == TM_OK) {
        status = tm_grow_array(
            (void **) &mesh->triangles,
            &mesh->triangle_capacity,
            sizeof(*mesh->triangles),
            mesh->triangle_count + 2
        );
    }
    if (status != TM_OK) {
        return status;
    }

    status = tm_append_point(mesh, point, -1, kind, &point_index);
    if (status != TM_OK) {
        return status;
    }

    if (!tm_make_triangle(&first, v0, v1, point_index, mesh->points) ||
        !tm_make_triangle(&second, v1, v2, point_index, mesh->points) ||
        !tm_make_triangle(&third, v2, v0, point_index, mesh->points)) {
        return TM_ERR_INTERNAL;
    }

    status = tm_copy_triangle_edge_state(&triangle, v0, v1, &first);
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&triangle, v1, v2, &second);
    }
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&triangle, v2, v0, &third);
    }
    if (status != TM_OK) {
        return status;
    }

    mesh->triangles[triangle_index] = first;
    mesh->triangles[second_index] = second;
    mesh->triangles[third_index] = third;
    mesh->triangle_count += 2;

    status = tm_link_triangles_across_edge(mesh, triangle_index, second_index, v1, point_index);
    if (status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, second_index, third_index, v2, point_index);
    }
    if (status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, third_index, triangle_index, v0, point_index);
    }
    if (status == TM_OK) {
        status = tm_redirect_neighbor_reference(
            mesh,
            triangle.nbr[tm_find_edge_in_triangle(&triangle, v1, v2)],
            v1,
            v2,
            triangle_index,
            second_index
        );
    }
    if (status == TM_OK) {
        status = tm_redirect_neighbor_reference(
            mesh,
            triangle.nbr[tm_find_edge_in_triangle(&triangle, v2, v0)],
            v2,
            v0,
            triangle_index,
            third_index
        );
    }
    if (status != TM_OK) {
        return status;
    }

    tm_assign_incident_triangle(mesh, v0, triangle_index);
    tm_assign_incident_triangle(mesh, v1, triangle_index);
    tm_assign_incident_triangle(mesh, v2, second_index);
    tm_assign_incident_triangle(mesh, point_index, triangle_index);
    seed_triangles[0] = triangle_index;
    seed_triangles[1] = second_index;
    seed_triangles[2] = third_index;
    status = tm_legalize_point_star(mesh, point_index, seed_triangles, 3);
    if (status != TM_OK) {
        return status;
    }
    if (out_point_index != NULL) {
        *out_point_index = point_index;
    }
    return TM_OK;
}

TMStatus tm_insert_point_on_edge(
    TMMesh *mesh,
    int triangle_index,
    int edge,
    const double point[2],
    TMVertexKind kind,
    int *out_point_index
)
{
    TMTriangle first;
    TMTriangle second;
    TMTriangle third;
    TMTriangle fourth;
    TMTriangle triangle;
    TMTriangle neighbor;
    int neighbor_index;
    int neighbor_edge = -1;
    int a;
    int b;
    int c;
    int d = -1;
    int point_index;
    int second_index;
    int fourth_index = -1;
    int seed_triangles[4];
    size_t seed_count = 0;
    unsigned char split_constrained;
    TMStatus status;

    if (mesh == NULL || point == NULL || triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count || edge < 0 || edge > 2) {
        tm_set_pslg_error_detail(
            "insert point on edge: invalid triangle=%d edge=%d",
            triangle_index,
            edge
        );
        return TM_ERR_INTERNAL;
    }

    tm_triangle_edge_vertices(&mesh->triangles[triangle_index], edge, &a, &b);
    c = mesh->triangles[triangle_index].v[edge];

    if (!tm_point_nearly_on_segment(&mesh->points[a], &mesh->points[b], point)) {
        tm_set_pslg_error_detail(
            "insert point on edge: point (%.12g, %.12g) is not on segment %d-%d",
            point[0],
            point[1],
            a,
            b
        );
        return TM_ERR_INTERNAL;
    }

    neighbor_index = mesh->triangles[triangle_index].nbr[edge];
    if (neighbor_index >= 0) {
        neighbor_edge = tm_find_edge_in_triangle(&mesh->triangles[neighbor_index], a, b);
        if (neighbor_edge < 0) {
            tm_set_pslg_error_detail(
                "insert point on edge: triangle %d edge %d has inconsistent neighbor %d",
                triangle_index,
                edge,
                neighbor_index
            );
            return TM_ERR_INVALID_MESH;
        }
        d = mesh->triangles[neighbor_index].v[neighbor_edge];
    }
    if (kind == TM_VERTEX_TRIANGLE_SPLIT) {
        if (tm_point_nearly_equals_vertex(&mesh->points[a], point)) {
            if (out_point_index != NULL) {
                *out_point_index = a;
            }
            return TM_OK;
        }
        if (tm_point_nearly_equals_vertex(&mesh->points[b], point)) {
            if (out_point_index != NULL) {
                *out_point_index = b;
            }
            return TM_OK;
        }
        if (tm_point_nearly_equals_vertex(&mesh->points[c], point)) {
            if (out_point_index != NULL) {
                *out_point_index = c;
            }
            return TM_OK;
        }
        if (d >= 0 && tm_point_nearly_equals_vertex(&mesh->points[d], point)) {
            if (out_point_index != NULL) {
                *out_point_index = d;
            }
            return TM_OK;
        }
    }
    triangle = mesh->triangles[triangle_index];
    if (neighbor_index >= 0) {
        neighbor = mesh->triangles[neighbor_index];
    }
    split_constrained = mesh->triangles[triangle_index].constrained[edge];
    second_index = (int) mesh->triangle_count;
    if (neighbor_index >= 0) {
        fourth_index = second_index + 1;
    }

    status = tm_grow_array((void **) &mesh->points, &mesh->point_capacity, sizeof(*mesh->points), mesh->point_count + 1);
    if (status == TM_OK) {
        status = tm_grow_array(
            (void **) &mesh->triangles,
            &mesh->triangle_capacity,
            sizeof(*mesh->triangles),
            mesh->triangle_count + ((neighbor_index >= 0) ? 2 : 1)
        );
    }
    if (status != TM_OK) {
        return status;
    }

    status = tm_append_point(mesh, point, -1, kind, &point_index);
    if (status != TM_OK) {
        return status;
    }

    if (!tm_make_triangle(&first, c, a, point_index, mesh->points) ||
        !tm_make_triangle(&second, c, point_index, b, mesh->points)) {
        tm_set_pslg_error_detail(
            "insert point on edge: splitting %d-%d produced degenerate primary triangles via point %d",
            a,
            b,
            point_index
        );
        return TM_ERR_INTERNAL;
    }

    if (neighbor_index >= 0 && status == TM_OK) {
        if (!tm_make_triangle(&third, d, b, point_index, mesh->points) ||
            !tm_make_triangle(&fourth, d, point_index, a, mesh->points)) {
            tm_set_pslg_error_detail(
                "insert point on edge: splitting %d-%d produced degenerate neighbor triangles via point %d",
                a,
                b,
                point_index
            );
            return TM_ERR_INTERNAL;
        }
    }
    if (status != TM_OK) {
        return status;
    }

    status = tm_copy_triangle_edge_state(&triangle, c, a, &first);
    if (status == TM_OK) {
        status = tm_copy_triangle_edge_state(&triangle, b, c, &second);
    }
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_copy_triangle_edge_state(&neighbor, d, b, &third);
    }
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_copy_triangle_edge_state(&neighbor, a, d, &fourth);
    }
    if (status != TM_OK) {
        return status;
    }

    mesh->triangles[triangle_index] = first;
    mesh->triangles[second_index] = second;
    mesh->triangle_count += 1;
    if (neighbor_index >= 0) {
        mesh->triangles[neighbor_index] = third;
        mesh->triangles[fourth_index] = fourth;
        mesh->triangle_count += 1;
    }

    status = tm_link_triangles_across_edge(mesh, triangle_index, second_index, c, point_index);
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, neighbor_index, fourth_index, d, point_index);
    }
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, triangle_index, fourth_index, a, point_index);
    }
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_link_triangles_across_edge(mesh, second_index, neighbor_index, point_index, b);
    }
    if (status == TM_OK) {
        status = tm_redirect_neighbor_reference(
            mesh,
            triangle.nbr[tm_find_edge_in_triangle(&triangle, b, c)],
            b,
            c,
            triangle_index,
            second_index
        );
    }
    if (neighbor_index >= 0 && status == TM_OK) {
        status = tm_redirect_neighbor_reference(
            mesh,
            neighbor.nbr[tm_find_edge_in_triangle(&neighbor, a, d)],
            a,
            d,
            neighbor_index,
            fourth_index
        );
    }
    if (status == TM_OK && split_constrained) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[triangle_index], a, point_index, 1);
    }
    if (status == TM_OK && split_constrained) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[second_index], point_index, b, 1);
    }
    if (neighbor_index >= 0 && status == TM_OK && split_constrained) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[neighbor_index], point_index, b, 1);
    }
    if (neighbor_index >= 0 && status == TM_OK && split_constrained) {
        status = tm_set_triangle_edge_constraint(&mesh->triangles[fourth_index], a, point_index, 1);
    }
    if (status != TM_OK) {
        tm_set_pslg_error_detail(
            "insert point on edge: local topology update failed after splitting %d-%d at point %d",
            a,
            b,
            point_index
        );
        return status;
    }

    tm_assign_incident_triangle(mesh, a, triangle_index);
    tm_assign_incident_triangle(mesh, b, second_index);
    tm_assign_incident_triangle(mesh, c, triangle_index);
    if (neighbor_index >= 0) {
        tm_assign_incident_triangle(mesh, d, neighbor_index);
    }
    tm_assign_incident_triangle(mesh, point_index, triangle_index);
    seed_triangles[seed_count++] = triangle_index;
    seed_triangles[seed_count++] = second_index;
    if (neighbor_index >= 0) {
        seed_triangles[seed_count++] = neighbor_index;
        seed_triangles[seed_count++] = fourth_index;
    }
    status = tm_legalize_point_star(mesh, point_index, seed_triangles, seed_count);
    if (status != TM_OK) {
        return status;
    }
    if (out_point_index != NULL) {
        *out_point_index = point_index;
    }
    return TM_OK;
}

void tm_free_points(TMPoint *points)
{
    free(points);
}

void tm_free_mesh(TMMesh *mesh)
{
    if (mesh == NULL) {
        return;
    }

    free(mesh->points);
    free(mesh->triangles);
    free(mesh->segments);
    mesh->points = NULL;
    mesh->triangles = NULL;
    mesh->segments = NULL;
    mesh->point_count = 0;
    mesh->point_capacity = 0;
    mesh->triangle_count = 0;
    mesh->triangle_capacity = 0;
    mesh->segment_count = 0;
    mesh->segment_capacity = 0;
}
