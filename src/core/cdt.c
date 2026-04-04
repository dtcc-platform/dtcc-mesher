#include "cdt.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "predicates.h"

static const double tm_pi = 3.14159265358979323846;
static const double tm_default_min_angle_deg = 20.0;

static int tm_find_live_segment_index(const TMMesh *mesh, int a, int b, size_t *out_segment_index);
static TMStatus tm_recover_segment(TMMesh *mesh, int a, int b, const TMBuildOptions *options);
static void tm_choose_recovery_split_coordinates(
    const TMMesh *mesh,
    int a,
    int b,
    double point[2],
    int *out_blocker_index,
    double *out_blocker_distance
);
static TMStatus tm_split_recovery_segment(
    TMMesh *mesh,
    int a,
    int b,
    const TMBuildOptions *options
);
static TMStatus tm_update_after_local_edit(TMMesh *mesh, const TMBuildOptions *options);
static TMStatus tm_restore_constrained_delaunay(TMMesh *mesh, const TMBuildOptions *options);

static void tm_verbose_log(const TMBuildOptions *options, const char *format, ...)
{
    va_list args;

    if (options == NULL || !options->verbose) {
        return;
    }

    va_start(args, format);
    vfprintf(stderr, format, args);
    fputc('\n', stderr);
    va_end(args);
}

static double tm_orient_xy(const double a[2], const double b[2], const double c[2])
{
    return orient2d((REAL *) a, (REAL *) b, (REAL *) c);
}

static int tm_refinement_enabled(const TMBuildOptions *options)
{
    return options != NULL && options->refine;
}

static int tm_use_offcenters(const TMBuildOptions *options)
{
    return options != NULL && options->use_offcenters;
}

static int tm_protect_acute_corners(const TMBuildOptions *options)
{
    return options == NULL || options->protect_acute_corners;
}

static TMAcuteProtectionMode tm_acute_mode(const TMBuildOptions *options)
{
    if (options != NULL) {
        return options->acute_mode;
    }

    return TM_ACUTE_MODE_SHELL;
}

static const char *tm_acute_mode_name(const TMBuildOptions *options)
{
    return tm_acute_mode(options) == TM_ACUTE_MODE_SIMPLE ? "simple" : "shell";
}

static double tm_quality_min_angle_deg(const TMBuildOptions *options)
{
    if (options != NULL && options->min_angle_deg > 0.0) {
        return options->min_angle_deg;
    }

    return tm_default_min_angle_deg;
}

static double tm_area_limit(const TMBuildOptions *options)
{
    if (options != NULL && options->max_area > 0.0) {
        return options->max_area;
    }

    return 0.0;
}

static double tm_protection_angle_deg(const TMBuildOptions *options)
{
    if (options != NULL && options->protect_angle_deg > 0.0) {
        return options->protect_angle_deg;
    }

    return tm_quality_min_angle_deg(options);
}

static size_t tm_refinement_step_limit(const TMBuildOptions *options, const TMMesh *mesh)
{
    if (options != NULL && options->max_refinement_steps != 0) {
        return options->max_refinement_steps;
    }

    return mesh->point_count * mesh->point_count * 8 + mesh->segment_count * 32 + 256;
}

static size_t tm_protection_level_limit(const TMBuildOptions *options)
{
    if (options != NULL && options->max_protection_levels != 0) {
        return options->max_protection_levels;
    }

    return 6;
}

static double tm_quality_beta(double min_angle_deg)
{
    return 1.0 / (2.0 * sin(min_angle_deg * tm_pi / 180.0));
}

static int tm_segments_cross_properly(const double a[2], const double b[2], const double c[2], const double d[2])
{
    double o1 = tm_orient_xy(a, b, c);
    double o2 = tm_orient_xy(a, b, d);
    double o3 = tm_orient_xy(c, d, a);
    double o4 = tm_orient_xy(c, d, b);

    return ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0)) &&
           ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0));
}

static int tm_point_nearly_on_segment_xy(const double a[2], const double b[2], const double point[2])
{
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double len_sq = dx * dx + dy * dy;
    double cross = dx * (point[1] - a[1]) - dy * (point[0] - a[0]);
    double dot = (point[0] - a[0]) * (point[0] - b[0]) +
                 (point[1] - a[1]) * (point[1] - b[1]);
    double tolerance;

    if (len_sq == 0.0) {
        return 0;
    }

    tolerance = 1e-12 * ((len_sq > 1.0) ? len_sq : 1.0);
    return fabs(cross) <= tolerance && dot <= tolerance;
}

static int tm_mesh_has_edge(const TMMesh *mesh, int a, int b, int *out_triangle, int *out_edge)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge = tm_find_edge_in_triangle(&mesh->triangles[tri_index], a, b);

        if (edge >= 0) {
            if (out_triangle != NULL) {
                *out_triangle = (int) tri_index;
            }
            if (out_edge != NULL) {
                *out_edge = edge;
            }
            return 1;
        }
    }

    return 0;
}

static TMStatus tm_mark_constraint_edge(TMMesh *mesh, int a, int b)
{
    size_t tri_index;
    int found = 0;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge = tm_find_edge_in_triangle(&mesh->triangles[tri_index], a, b);

        if (edge >= 0) {
            mesh->triangles[tri_index].constrained[edge] = 1;
            found = 1;
        }
    }

    return found ? TM_OK : TM_ERR_INVALID_MESH;
}

static TMStatus tm_mark_all_constraints(TMMesh *mesh)
{
    size_t tri_index;
    size_t segment_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            mesh->triangles[tri_index].constrained[edge] = 0;
        }
    }

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        TMStatus status;

        if (!mesh->segments[segment_index].live) {
            continue;
        }
        status = tm_mark_constraint_edge(
            mesh,
            mesh->segments[segment_index].v[0],
            mesh->segments[segment_index].v[1]
        );
        if (status != TM_OK) {
            tm_set_pslg_error_detail(
                "live segment %d-%d is not present as a mesh edge after constraint recovery",
                mesh->segments[segment_index].v[0],
                mesh->segments[segment_index].v[1]
            );
            return status;
        }
    }

    return TM_OK;
}

static TMStatus tm_copy_segments_to_mesh(const TMPSLG *pslg, TMMesh *mesh)
{
    size_t i;

    mesh->segments = (TMSegment *) calloc(pslg->segment_count, sizeof(*mesh->segments));
    if (mesh->segments == NULL && pslg->segment_count != 0) {
        return TM_ERR_ALLOC;
    }

    mesh->segment_count = pslg->segment_count;
    mesh->segment_capacity = pslg->segment_count;

    for (i = 0; i < pslg->segment_count; ++i) {
        mesh->segments[i] = pslg->segments[i];
        mesh->segments[i].live = 1;
        mesh->segments[i].is_protected = 0;
        mesh->segments[i].protected_apex = -1;
    }

    return TM_OK;
}

static double tm_distance_squared(const double a[2], const double b[2])
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];

    return dx * dx + dy * dy;
}

static double tm_segment_length_squared(const TMMesh *mesh, int a, int b)
{
    return tm_distance_squared(mesh->points[a].xy, mesh->points[b].xy);
}

static TMStatus tm_clone_pslg(const TMPSLG *source, TMPSLG *target)
{
    memset(target, 0, sizeof(*target));

    if (source->point_count != 0) {
        target->points = (TMPoint *) calloc(source->point_count, sizeof(*target->points));
        if (target->points == NULL) {
            return TM_ERR_ALLOC;
        }
        memcpy(target->points, source->points, source->point_count * sizeof(*target->points));
        target->point_count = source->point_count;
    }

    if (source->segment_count != 0) {
        target->segments = (TMSegment *) calloc(source->segment_count, sizeof(*target->segments));
        if (target->segments == NULL) {
            tm_free_pslg(target);
            return TM_ERR_ALLOC;
        }
        memcpy(target->segments, source->segments, source->segment_count * sizeof(*target->segments));
        target->segment_count = source->segment_count;
    }

    if (source->hole_count != 0) {
        target->holes = (double (*)[2]) calloc(source->hole_count, sizeof(*target->holes));
        if (target->holes == NULL) {
            tm_free_pslg(target);
            return TM_ERR_ALLOC;
        }
        memcpy(target->holes, source->holes, source->hole_count * sizeof(*target->holes));
        target->hole_count = source->hole_count;
    }

    return TM_OK;
}

static TMStatus tm_split_pslg_segment_at_point(TMPSLG *pslg, size_t segment_index, const double point[2])
{
    TMSegment segment = pslg->segments[segment_index];
    TMPoint *points;
    TMSegment *segments;
    int point_index = (int) pslg->point_count;

    points = (TMPoint *) realloc(pslg->points, (pslg->point_count + 1) * sizeof(*points));
    if (points == NULL) {
        return TM_ERR_ALLOC;
    }
    pslg->points = points;
    pslg->points[point_index].xy[0] = point[0];
    pslg->points[point_index].xy[1] = point[1];
    pslg->points[point_index].original_index = point_index;
    pslg->points[point_index].kind = TM_VERTEX_SEGMENT_SPLIT;
    pslg->points[point_index].incident_triangle = -1;
    pslg->points[point_index].protection_apex = -1;
    pslg->points[point_index].protection_side = TM_PROTECTION_SIDE_NONE;
    pslg->points[point_index].protection_level = 0;
    pslg->point_count += 1;

    segments = (TMSegment *) realloc(pslg->segments, (pslg->segment_count + 1) * sizeof(*segments));
    if (segments == NULL) {
        return TM_ERR_ALLOC;
    }
    pslg->segments = segments;
    pslg->segments[segment_index].v[1] = point_index;
    pslg->segments[pslg->segment_count] = segment;
    pslg->segments[pslg->segment_count].v[0] = point_index;
    pslg->segment_count += 1;
    return TM_OK;
}

static TMStatus tm_append_pslg_point(TMPSLG *pslg, const double point[2], TMVertexKind kind)
{
    TMPoint *points = (TMPoint *) realloc(pslg->points, (pslg->point_count + 1) * sizeof(*points));
    size_t point_index = pslg->point_count;

    if (points == NULL) {
        return TM_ERR_ALLOC;
    }

    pslg->points = points;
    pslg->points[point_index].xy[0] = point[0];
    pslg->points[point_index].xy[1] = point[1];
    pslg->points[point_index].original_index = (int) point_index;
    pslg->points[point_index].kind = kind;
    pslg->points[point_index].incident_triangle = -1;
    pslg->points[point_index].protection_apex = -1;
    pslg->points[point_index].protection_side = TM_PROTECTION_SIDE_NONE;
    pslg->points[point_index].protection_level = 0;
    pslg->point_count += 1;
    return TM_OK;
}

static double tm_segment_length_limit(const TMBuildOptions *options)
{
    if (options != NULL && options->max_edge_length > 0.0) {
        return options->max_edge_length;
    }

    double max_area = tm_area_limit(options);

    if (max_area <= 0.0) {
        return 0.0;
    }

    return sqrt(2.0 * max_area);
}

static double tm_point_segment_distance_squared(const double point[2], const double a[2], const double b[2])
{
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double length_sq = dx * dx + dy * dy;

    if (length_sq <= 0.0) {
        return tm_distance_squared(point, a);
    }

    {
        double t = ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / length_sq;
        double projection[2];

        if (t < 0.0) {
            t = 0.0;
        } else if (t > 1.0) {
            t = 1.0;
        }

        projection[0] = a[0] + t * dx;
        projection[1] = a[1] + t * dy;
        return tm_distance_squared(point, projection);
    }
}

static int tm_find_near_segment_free_point(
    const TMPSLG *pslg,
    size_t segment_index,
    const size_t *degree,
    int *out_point_index,
    double out_projection[2],
    double *out_distance
)
{
    const TMSegment *segment = &pslg->segments[segment_index];
    const double *a = pslg->points[segment->v[0]].xy;
    const double *b = pslg->points[segment->v[1]].xy;
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double segment_length_sq = dx * dx + dy * dy;
    double best_distance_sq = INFINITY;
    double best_t = 0.0;
    int best_point_index = -1;
    size_t point_index;

    if (segment_length_sq <= 0.0) {
        return 0;
    }

    for (point_index = 0; point_index < pslg->point_count; ++point_index) {
        const double *candidate = pslg->points[point_index].xy;
        double t;
        double projection[2];
        double distance_sq;

        if (degree[point_index] != 0 || (int) point_index == segment->v[0] || (int) point_index == segment->v[1]) {
            continue;
        }

        t = ((candidate[0] - a[0]) * dx + (candidate[1] - a[1]) * dy) / segment_length_sq;
        if (t <= 1e-6 || t >= 1.0 - 1e-6) {
            continue;
        }

        projection[0] = a[0] + t * dx;
        projection[1] = a[1] + t * dy;
        distance_sq = tm_distance_squared(candidate, projection);
        if (distance_sq < best_distance_sq) {
            best_distance_sq = distance_sq;
            best_t = t;
            best_point_index = (int) point_index;
        }
    }

    if (best_point_index < 0 || best_distance_sq > segment_length_sq * 1e-4) {
        return 0;
    }

    out_projection[0] = a[0] + best_t * dx;
    out_projection[1] = a[1] + best_t * dy;
    if (out_point_index != NULL) {
        *out_point_index = best_point_index;
    }
    if (out_distance != NULL) {
        *out_distance = sqrt(best_distance_sq);
    }
    return 1;
}

static TMStatus tm_pre_split_near_segment_free_points(TMPSLG *pslg, const TMBuildOptions *options)
{
    for (;;) {
        size_t *degree;
        size_t segment_index;
        int changed = 0;

        degree = (size_t *) calloc(pslg->point_count, sizeof(*degree));
        if (degree == NULL) {
            return TM_ERR_ALLOC;
        }

        for (segment_index = 0; segment_index < pslg->segment_count; ++segment_index) {
            degree[pslg->segments[segment_index].v[0]] += 1;
            degree[pslg->segments[segment_index].v[1]] += 1;
        }

        for (segment_index = 0; segment_index < pslg->segment_count; ++segment_index) {
            int blocker_point = -1;
            double projection[2];
            double blocker_distance = 0.0;
            TMStatus status;

            if (!tm_find_near_segment_free_point(pslg, segment_index, degree, &blocker_point, projection, &blocker_distance)) {
                continue;
            }

            tm_verbose_log(
                options,
                "pre-split segment %d-%d at projection of free point %d (distance %.6g)",
                pslg->segments[segment_index].v[0],
                pslg->segments[segment_index].v[1],
                blocker_point,
                blocker_distance
            );
            status = tm_split_pslg_segment_at_point(pslg, segment_index, projection);
            free(degree);
            if (status != TM_OK) {
                return status;
            }
            changed = 1;
            break;
        }

        if (!changed) {
            free(degree);
            return TM_OK;
        }
    }
}

static TMStatus tm_pre_split_long_segments(TMPSLG *pslg, const TMBuildOptions *options)
{
    double max_length = tm_segment_length_limit(options);

    if (max_length <= 0.0) {
        return TM_OK;
    }

    for (;;) {
        size_t segment_index;
        int changed = 0;

        for (segment_index = 0; segment_index < pslg->segment_count; ++segment_index) {
            const TMSegment *segment = &pslg->segments[segment_index];
            const double *a = pslg->points[segment->v[0]].xy;
            const double *b = pslg->points[segment->v[1]].xy;
            double length = sqrt(tm_distance_squared(a, b));
            double midpoint[2];
            TMStatus status;

            if (length <= max_length * (1.0 + 1e-12)) {
                continue;
            }

            midpoint[0] = 0.5 * (a[0] + b[0]);
            midpoint[1] = 0.5 * (a[1] + b[1]);
            tm_verbose_log(
                options,
                "pre-split long segment %d-%d (length %.6g > %.6g)",
                segment->v[0],
                segment->v[1],
                length,
                max_length
            );
            status = tm_split_pslg_segment_at_point(pslg, segment_index, midpoint);
            if (status != TM_OK) {
                return status;
            }
            changed = 1;
            break;
        }

        if (!changed) {
            return TM_OK;
        }
    }
}

static int tm_point_in_domain(const TMPSLG *pslg, const double point[2]);

static TMStatus tm_seed_interior_pslg_points(TMPSLG *pslg, const TMBuildOptions *options)
{
    double spacing = tm_segment_length_limit(options);
    double vertical_step;
    double boundary_clearance;
    double minx;
    double miny;
    double maxx;
    double maxy;
    double clearance_sq;
    size_t point_index;
    size_t added = 0;
    int row = 0;
    double y;

    if (spacing <= 0.0 || pslg->point_count == 0) {
        return TM_OK;
    }

    minx = pslg->points[0].xy[0];
    miny = pslg->points[0].xy[1];
    maxx = pslg->points[0].xy[0];
    maxy = pslg->points[0].xy[1];
    for (point_index = 1; point_index < pslg->point_count; ++point_index) {
        const double *point = pslg->points[point_index].xy;

        if (point[0] < minx) {
            minx = point[0];
        }
        if (point[1] < miny) {
            miny = point[1];
        }
        if (point[0] > maxx) {
            maxx = point[0];
        }
        if (point[1] > maxy) {
            maxy = point[1];
        }
    }

    vertical_step = spacing * 0.8660254037844386;
    boundary_clearance = 0.5 * vertical_step;
    clearance_sq = boundary_clearance * boundary_clearance;

    y = miny + 0.5 * vertical_step;
    while (y < maxy) {
        double x_offset = (row % 2) ? 0.5 * spacing : 0.0;
        double x = minx + 0.5 * spacing + x_offset;

        while (x < maxx) {
            double candidate[2];
            double min_distance_sq = INFINITY;
            size_t segment_index;

            candidate[0] = x;
            candidate[1] = y;
            if (tm_point_in_domain(pslg, candidate)) {
                for (segment_index = 0; segment_index < pslg->segment_count; ++segment_index) {
                    const double *a = pslg->points[pslg->segments[segment_index].v[0]].xy;
                    const double *b = pslg->points[pslg->segments[segment_index].v[1]].xy;
                    double distance_sq = tm_point_segment_distance_squared(candidate, a, b);

                    if (distance_sq < min_distance_sq) {
                        min_distance_sq = distance_sq;
                    }
                    if (min_distance_sq < clearance_sq) {
                        break;
                    }
                }

                if (min_distance_sq >= clearance_sq) {
                    TMStatus status = tm_append_pslg_point(pslg, candidate, TM_VERTEX_INPUT);

                    if (status != TM_OK) {
                        return status;
                    }
                    added += 1;
                }
            }

            x += spacing;
        }

        y += vertical_step;
        row += 1;
    }

    if (added > 0) {
        tm_verbose_log(options, "seeded %zu interior size points", added);
    }
    return TM_OK;
}

static int tm_point_encroaches_segment(const double point[2], const double a[2], const double b[2], double segment_length_sq)
{
    double dot = (point[0] - a[0]) * (point[0] - b[0]) + (point[1] - a[1]) * (point[1] - b[1]);
    double tolerance = 1e-12 * ((segment_length_sq > 1.0) ? segment_length_sq : 1.0);

    return dot < -tolerance;
}

static TMStatus tm_grow_segments(TMMesh *mesh, size_t min_capacity)
{
    TMSegment *segments;
    size_t new_capacity;

    if (mesh->segment_capacity >= min_capacity) {
        return TM_OK;
    }

    new_capacity = (mesh->segment_capacity == 0) ? 8 : mesh->segment_capacity;
    while (new_capacity < min_capacity) {
        if (new_capacity > ((size_t) -1) / 2) {
            return TM_ERR_ALLOC;
        }
        new_capacity *= 2;
    }

    segments = (TMSegment *) realloc(mesh->segments, new_capacity * sizeof(*segments));
    if (segments == NULL) {
        return TM_ERR_ALLOC;
    }

    mesh->segments = segments;
    mesh->segment_capacity = new_capacity;
    return TM_OK;
}

static TMStatus tm_append_live_segment(TMMesh *mesh, int a, int b, int original_index)
{
    TMStatus status;

    status = tm_grow_segments(mesh, mesh->segment_count + 1);
    if (status != TM_OK) {
        return status;
    }

    mesh->segments[mesh->segment_count].v[0] = a;
    mesh->segments[mesh->segment_count].v[1] = b;
    mesh->segments[mesh->segment_count].original_index = original_index;
    mesh->segments[mesh->segment_count].live = 1;
    mesh->segments[mesh->segment_count].is_protected = 0;
    mesh->segments[mesh->segment_count].protected_apex = -1;
    mesh->segment_count += 1;
    return TM_OK;
}

static int tm_find_live_incident_segments(const TMMesh *mesh, int apex, size_t out_segments[2])
{
    size_t segment_index;
    size_t count = 0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        if (!mesh->segments[segment_index].live) {
            continue;
        }
        if (mesh->segments[segment_index].v[0] == apex || mesh->segments[segment_index].v[1] == apex) {
            if (count >= 2) {
                return 0;
            }
            out_segments[count] = segment_index;
            count += 1;
        }
    }

    return count == 2;
}

static int tm_segment_other_endpoint(const TMSegment *segment, int vertex)
{
    return (segment->v[0] == vertex) ? segment->v[1] : segment->v[0];
}

static double tm_corner_angle_deg(const TMMesh *mesh, int apex, int left, int right)
{
    double ux = mesh->points[left].xy[0] - mesh->points[apex].xy[0];
    double uy = mesh->points[left].xy[1] - mesh->points[apex].xy[1];
    double vx = mesh->points[right].xy[0] - mesh->points[apex].xy[0];
    double vy = mesh->points[right].xy[1] - mesh->points[apex].xy[1];
    double unorm = sqrt(ux * ux + uy * uy);
    double vnorm = sqrt(vx * vx + vy * vy);
    double cosine;

    if (unorm <= 0.0 || vnorm <= 0.0) {
        return 180.0;
    }

    cosine = (ux * vx + uy * vy) / (unorm * vnorm);
    if (cosine < -1.0) {
        cosine = -1.0;
    }
    if (cosine > 1.0) {
        cosine = 1.0;
    }

    return acos(cosine) * 180.0 / tm_pi;
}

static int tm_local_triangle_is_exempt_from_quality(const TMMesh *mesh, int triangle_index)
{
    return tm_triangle_is_quality_exempt(mesh, triangle_index);
}

static int tm_find_live_segment_index(const TMMesh *mesh, int a, int b, size_t *out_segment_index)
{
    size_t segment_index;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        int seg_a;
        int seg_b;
        int seg_lo;
        int seg_hi;

        if (!mesh->segments[segment_index].live) {
            continue;
        }

        seg_a = mesh->segments[segment_index].v[0];
        seg_b = mesh->segments[segment_index].v[1];
        seg_lo = (seg_a < seg_b) ? seg_a : seg_b;
        seg_hi = (seg_a < seg_b) ? seg_b : seg_a;
        if (seg_lo == lo && seg_hi == hi) {
            if (out_segment_index != NULL) {
                *out_segment_index = segment_index;
            }
            return 1;
        }
    }

    return 0;
}

static TMStatus tm_split_segment_registry(TMMesh *mesh, size_t segment_index, int point_index)
{
    TMSegment segment = mesh->segments[segment_index];
    TMStatus status;
    size_t first_new_index;
    size_t second_new_index;

    mesh->segments[segment_index].live = 0;

    first_new_index = mesh->segment_count;
    status = tm_append_live_segment(mesh, segment.v[0], point_index, segment.original_index);
    if (status == TM_OK) {
        mesh->segments[first_new_index].is_protected = segment.is_protected;
        mesh->segments[first_new_index].protected_apex = segment.protected_apex;
        second_new_index = mesh->segment_count;
        status = tm_append_live_segment(mesh, point_index, segment.v[1], segment.original_index);
        if (status == TM_OK) {
            mesh->segments[second_new_index].is_protected = segment.is_protected;
            mesh->segments[second_new_index].protected_apex = segment.protected_apex;
        }
    }
    return status;
}

static void tm_copy_segment_flags(TMSegment *destination, const TMSegment *source)
{
    destination->is_protected = source->is_protected;
    destination->protected_apex = source->protected_apex;
}

static TMStatus tm_copy_live_segments_with_split(
    const TMMesh *source,
    TMMesh *target,
    int a,
    int b,
    int point_index
)
{
    size_t segment_index;

    for (segment_index = 0; segment_index < source->segment_count; ++segment_index) {
        const TMSegment *segment = &source->segments[segment_index];
        int matches_target =
            segment->live &&
            ((segment->v[0] == a && segment->v[1] == b) || (segment->v[0] == b && segment->v[1] == a));

        if (!segment->live) {
            continue;
        }

        if (matches_target) {
            size_t first_index = target->segment_count;
            TMStatus status = tm_append_live_segment(target, a, point_index, segment->original_index);
            if (status != TM_OK) {
                return status;
            }
            tm_copy_segment_flags(&target->segments[first_index], segment);

            first_index = target->segment_count;
            status = tm_append_live_segment(target, point_index, b, segment->original_index);
            if (status != TM_OK) {
                return status;
            }
            tm_copy_segment_flags(&target->segments[first_index], segment);
            continue;
        }

        {
            size_t new_index = target->segment_count;
            TMStatus status = tm_append_live_segment(target, segment->v[0], segment->v[1], segment->original_index);
            if (status != TM_OK) {
                return status;
            }
            tm_copy_segment_flags(&target->segments[new_index], segment);
        }
    }

    return TM_OK;
}

static TMStatus tm_recover_all_live_segments(TMMesh *mesh, const TMBuildOptions *options)
{
    size_t segment_index;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        TMSegment segment = mesh->segments[segment_index];
        TMStatus status;

        if (!segment.live) {
            continue;
        }

        status = tm_recover_segment(mesh, segment.v[0], segment.v[1], options);
        if (status != TM_OK) {
            return status;
        }
    }

    return tm_restore_constrained_delaunay(mesh, options);
}

static double tm_triangle_radius_edge_ratio(const TMMesh *mesh, int triangle_index)
{
    const TMTriangle *triangle = &mesh->triangles[triangle_index];
    double len01 = sqrt(tm_distance_squared(mesh->points[triangle->v[0]].xy, mesh->points[triangle->v[1]].xy));
    double len12 = sqrt(tm_distance_squared(mesh->points[triangle->v[1]].xy, mesh->points[triangle->v[2]].xy));
    double len20 = sqrt(tm_distance_squared(mesh->points[triangle->v[2]].xy, mesh->points[triangle->v[0]].xy));
    double shortest = len01;
    double area = fabs(orient2d(
        (REAL *) mesh->points[triangle->v[0]].xy,
        (REAL *) mesh->points[triangle->v[1]].xy,
        (REAL *) mesh->points[triangle->v[2]].xy
    )) * 0.5;

    if (len12 < shortest) {
        shortest = len12;
    }
    if (len20 < shortest) {
        shortest = len20;
    }
    if (shortest <= 0.0 || area <= 0.0) {
        return INFINITY;
    }

    return (len01 * len12 * len20) / (4.0 * area * shortest);
}

static double tm_triangle_area(const TMMesh *mesh, int triangle_index)
{
    const TMTriangle *triangle = &mesh->triangles[triangle_index];

    return fabs(orient2d(
        (REAL *) mesh->points[triangle->v[0]].xy,
        (REAL *) mesh->points[triangle->v[1]].xy,
        (REAL *) mesh->points[triangle->v[2]].xy
    )) * 0.5;
}

static TMStatus tm_triangle_shortest_edge(
    const TMMesh *mesh,
    int triangle_index,
    int *out_edge,
    int *out_a,
    int *out_b,
    int *out_apex,
    double *out_length_sq
)
{
    const TMTriangle *triangle = &mesh->triangles[triangle_index];
    int edge;
    int best_edge = 0;
    double best_length_sq = -1.0;

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double length_sq;

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        length_sq = tm_segment_length_squared(mesh, a, b);
        if (best_length_sq < 0.0 || length_sq < best_length_sq) {
            best_length_sq = length_sq;
            best_edge = edge;
            if (out_a != NULL) {
                *out_a = a;
            }
            if (out_b != NULL) {
                *out_b = b;
            }
            if (out_apex != NULL) {
                *out_apex = triangle->v[edge];
            }
        }
    }

    if (best_length_sq <= 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    if (out_edge != NULL) {
        *out_edge = best_edge;
    }
    if (out_length_sq != NULL) {
        *out_length_sq = best_length_sq;
    }
    return TM_OK;
}

static TMStatus tm_triangle_circumcenter(const TMMesh *mesh, int triangle_index, double out_point[2])
{
    const TMTriangle *triangle = &mesh->triangles[triangle_index];
    const double *a = mesh->points[triangle->v[0]].xy;
    const double *b = mesh->points[triangle->v[1]].xy;
    const double *c = mesh->points[triangle->v[2]].xy;
    double aa = a[0] * a[0] + a[1] * a[1];
    double bb = b[0] * b[0] + b[1] * b[1];
    double cc = c[0] * c[0] + c[1] * c[1];
    double det = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));

    if (det == 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    out_point[0] = (aa * (b[1] - c[1]) + bb * (c[1] - a[1]) + cc * (a[1] - b[1])) / det;
    out_point[1] = (aa * (c[0] - b[0]) + bb * (a[0] - c[0]) + cc * (b[0] - a[0])) / det;
    return TM_OK;
}

static int tm_triangle_offcenter_candidate(
    const TMMesh *mesh,
    int triangle_index,
    double min_angle_deg,
    const double circumcenter[2],
    double out_point[2]
)
{
    int a;
    int b;
    int apex;
    double length_sq;
    double length;
    double midpoint[2];
    double normal[2];
    double to_apex[2];
    double off_distance;
    double off_distance_sq;
    double circumcenter_distance_sq;
    double tangent;

    if (tm_triangle_shortest_edge(mesh, triangle_index, NULL, &a, &b, &apex, &length_sq) != TM_OK) {
        return 0;
    }

    length = sqrt(length_sq);
    tangent = tan(min_angle_deg * tm_pi / 360.0);
    if (length <= 0.0 || tangent <= 0.0) {
        return 0;
    }

    midpoint[0] = 0.5 * (mesh->points[a].xy[0] + mesh->points[b].xy[0]);
    midpoint[1] = 0.5 * (mesh->points[a].xy[1] + mesh->points[b].xy[1]);
    normal[0] = -(mesh->points[b].xy[1] - mesh->points[a].xy[1]) / length;
    normal[1] = (mesh->points[b].xy[0] - mesh->points[a].xy[0]) / length;
    to_apex[0] = mesh->points[apex].xy[0] - midpoint[0];
    to_apex[1] = mesh->points[apex].xy[1] - midpoint[1];
    if (normal[0] * to_apex[0] + normal[1] * to_apex[1] < 0.0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
    }

    off_distance = length / (2.0 * tangent);
    off_distance_sq = off_distance * off_distance;
    circumcenter_distance_sq = tm_distance_squared(midpoint, circumcenter);
    if (off_distance_sq > circumcenter_distance_sq * (1.0 + 1e-12)) {
        return 0;
    }

    out_point[0] = midpoint[0] + normal[0] * off_distance;
    out_point[1] = midpoint[1] + normal[1] * off_distance;
    return 1;
}

static int tm_point_inside_or_on_triangle(const TMMesh *mesh, int triangle_index, const double point[2], int *out_edge)
{
    const TMTriangle *triangle = &mesh->triangles[triangle_index];
    int edge;
    int on_edge = -1;

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double side;
        int near_edge;

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        side = tm_orient_xy(mesh->points[a].xy, mesh->points[b].xy, point);
        near_edge = tm_point_nearly_on_segment_xy(mesh->points[a].xy, mesh->points[b].xy, point);
        if (side < 0.0) {
            if (near_edge) {
                on_edge = edge;
                continue;
            }
            return 0;
        }
        if (side == 0.0 || near_edge) {
            on_edge = edge;
        }
    }

    if (out_edge != NULL) {
        *out_edge = on_edge;
    }
    return 1;
}

static int tm_find_encroached_segment(const TMMesh *mesh, size_t *out_segment_index, int *out_point_index)
{
    size_t segment_index;
    int found = 0;
    double best_length_sq = 0.0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &mesh->segments[segment_index];
        const double *a;
        const double *b;
        double length_sq;
        size_t point_index;

        if (!segment->live) {
            continue;
        }
        if (segment->is_protected) {
            continue;
        }

        a = mesh->points[segment->v[0]].xy;
        b = mesh->points[segment->v[1]].xy;
        length_sq = tm_segment_length_squared(mesh, segment->v[0], segment->v[1]);

        for (point_index = 0; point_index < mesh->point_count; ++point_index) {
            if ((int) point_index == segment->v[0] || (int) point_index == segment->v[1]) {
                continue;
            }
            if (mesh->points[point_index].incident_triangle < 0) {
                continue;
            }

            if (tm_point_encroaches_segment(mesh->points[point_index].xy, a, b, length_sq)) {
                if (!found ||
                    length_sq < best_length_sq ||
                    (length_sq == best_length_sq && segment_index < *out_segment_index)) {
                    found = 1;
                    best_length_sq = length_sq;
                    *out_segment_index = segment_index;
                    if (out_point_index != NULL) {
                        *out_point_index = (int) point_index;
                    }
                }
                break;
            }
        }
    }

    return found;
}

static int tm_find_encroached_segment_for_point(const TMMesh *mesh, const double point[2], size_t *out_segment_index)
{
    size_t segment_index;
    int found = 0;
    double best_length_sq = 0.0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &mesh->segments[segment_index];
        double length_sq;

        if (!segment->live) {
            continue;
        }
        if (segment->is_protected) {
            continue;
        }

        length_sq = tm_segment_length_squared(mesh, segment->v[0], segment->v[1]);
        if (tm_point_encroaches_segment(
                point,
                mesh->points[segment->v[0]].xy,
                mesh->points[segment->v[1]].xy,
                length_sq
            )) {
            if (!found ||
                length_sq < best_length_sq ||
                (length_sq == best_length_sq && segment_index < *out_segment_index)) {
                found = 1;
                best_length_sq = length_sq;
                *out_segment_index = segment_index;
            }
        }
    }

    return found;
}

static int tm_find_bad_triangle(
    const TMMesh *mesh,
    double beta,
    double max_area,
    int shortest_edge_first,
    size_t *out_triangle_index,
    double *out_ratio,
    double *out_area,
    int *out_quality_violation
)
{
    size_t tri_index;
    int found_quality = 0;
    int found_area = 0;
    double best_ratio = beta;
    double best_shortest_edge_sq = 0.0;
    double best_area = max_area;
    size_t best_area_index = 0;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        double ratio = tm_triangle_radius_edge_ratio(mesh, (int) tri_index);
        int quality_exempt = tm_local_triangle_is_exempt_from_quality(mesh, (int) tri_index);

        if (!quality_exempt && ratio > beta * (1.0 + 1e-12)) {
            if (shortest_edge_first) {
                double shortest_edge_sq;

                if (tm_triangle_shortest_edge(mesh, (int) tri_index, NULL, NULL, NULL, NULL, &shortest_edge_sq) != TM_OK) {
                    continue;
                }

                if (!found_quality ||
                    shortest_edge_sq < best_shortest_edge_sq ||
                    (shortest_edge_sq == best_shortest_edge_sq && ratio > best_ratio)) {
                    found_quality = 1;
                    best_shortest_edge_sq = shortest_edge_sq;
                    best_ratio = ratio;
                    *out_triangle_index = tri_index;
                }
            } else if (!found_quality || ratio > best_ratio) {
                found_quality = 1;
                best_ratio = ratio;
                *out_triangle_index = tri_index;
            }
        }

        if (max_area > 0.0 && !quality_exempt) {
            double area = tm_triangle_area(mesh, (int) tri_index);

            if (area > max_area * (1.0 + 1e-12) &&
                (!found_area || area > best_area || (area == best_area && tri_index < best_area_index))) {
                found_area = 1;
                best_area = area;
                best_area_index = tri_index;
            }
        }
    }

    if (found_quality) {
        if (out_ratio != NULL) {
            *out_ratio = best_ratio;
        }
        if (out_area != NULL) {
            *out_area = tm_triangle_area(mesh, (int) *out_triangle_index);
        }
        if (out_quality_violation != NULL) {
            *out_quality_violation = 1;
        }
        return 1;
    }

    if (found_area) {
        *out_triangle_index = best_area_index;
        if (out_ratio != NULL) {
            *out_ratio = tm_triangle_radius_edge_ratio(mesh, (int) best_area_index);
        }
        if (out_area != NULL) {
            *out_area = best_area;
        }
        if (out_quality_violation != NULL) {
            *out_quality_violation = 0;
        }
        return 1;
    }

    if (out_ratio != NULL) {
        *out_ratio = best_ratio;
    }
    return 0;
}

static size_t tm_count_live_segments(const TMMesh *mesh)
{
    size_t segment_index;
    size_t live_count = 0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        if (mesh->segments[segment_index].live) {
            live_count += 1;
        }
    }

    return live_count;
}

static TMStatus tm_find_crossing_edge(
    const TMMesh *mesh,
    int seg_a,
    int seg_b,
    int *out_triangle,
    int *out_edge
)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            int a;
            int b;
            int neighbor = mesh->triangles[tri_index].nbr[edge];

            if (neighbor < 0 || (size_t) neighbor <= tri_index || mesh->triangles[tri_index].constrained[edge]) {
                continue;
            }

            tm_triangle_edge_vertices(&mesh->triangles[tri_index], edge, &a, &b);
            if (a == seg_a || a == seg_b || b == seg_a || b == seg_b) {
                continue;
            }

            if (tm_segments_cross_properly(
                    mesh->points[seg_a].xy,
                    mesh->points[seg_b].xy,
                    mesh->points[a].xy,
                    mesh->points[b].xy
                )) {
                *out_triangle = (int) tri_index;
                *out_edge = edge;
                return TM_OK;
            }
        }
    }

    return TM_ERR_INVALID_MESH;
}

static void tm_choose_recovery_split_coordinates(
    const TMMesh *mesh,
    int a,
    int b,
    double point[2],
    int *out_blocker_index,
    double *out_blocker_distance
)
{
    const double *pa = mesh->points[a].xy;
    const double *pb = mesh->points[b].xy;
    double dx = pb[0] - pa[0];
    double dy = pb[1] - pa[1];
    double segment_length_sq = dx * dx + dy * dy;
    double best_distance_sq = INFINITY;
    double best_t = 0.5;
    int best_index = -1;
    size_t point_index;

    point[0] = 0.5 * (pa[0] + pb[0]);
    point[1] = 0.5 * (pa[1] + pb[1]);

    if (segment_length_sq <= 0.0) {
        if (out_blocker_index != NULL) {
            *out_blocker_index = -1;
        }
        if (out_blocker_distance != NULL) {
            *out_blocker_distance = 0.0;
        }
        return;
    }

    for (point_index = 0; point_index < mesh->point_count; ++point_index) {
        const double *candidate = mesh->points[point_index].xy;
        double t;
        double projection[2];
        double distance_sq;

        if ((int) point_index == a || (int) point_index == b || mesh->points[point_index].kind == TM_VERTEX_SUPER) {
            continue;
        }

        t = ((candidate[0] - pa[0]) * dx + (candidate[1] - pa[1]) * dy) / segment_length_sq;
        if (t <= 1e-6 || t >= 1.0 - 1e-6) {
            continue;
        }

        projection[0] = pa[0] + t * dx;
        projection[1] = pa[1] + t * dy;
        distance_sq = tm_distance_squared(candidate, projection);
        if (distance_sq < best_distance_sq) {
            best_distance_sq = distance_sq;
            best_t = t;
            best_index = (int) point_index;
        }
    }

    if (best_index >= 0 && best_distance_sq <= segment_length_sq * 1e-6) {
        point[0] = pa[0] + best_t * dx;
        point[1] = pa[1] + best_t * dy;
    } else {
        best_index = -1;
        best_distance_sq = tm_distance_squared(pa, point);
    }

    if (out_blocker_index != NULL) {
        *out_blocker_index = best_index;
    }
    if (out_blocker_distance != NULL) {
        *out_blocker_distance = sqrt(best_distance_sq);
    }
}

static TMStatus tm_split_recovery_segment(
    TMMesh *mesh,
    int a,
    int b,
    const TMBuildOptions *options
)
{
    TMPoint *input_points = NULL;
    TMMesh rebuilt;
    double split_point[2];
    int blocker_index = -1;
    double blocker_distance = 0.0;
    int point_index;
    TMStatus status;

    tm_verbose_log(options, "recover segment %d-%d stalled; split and retry", a, b);
    tm_choose_recovery_split_coordinates(mesh, a, b, split_point, &blocker_index, &blocker_distance);
    if (blocker_index >= 0) {
        tm_verbose_log(
            options,
            "recover segment %d-%d: rebuild with split at projection of point %d (distance %.6g)",
            a,
            b,
            blocker_index,
            blocker_distance
        );
    } else {
        tm_verbose_log(options, "recover segment %d-%d: rebuild with midpoint split", a, b);
    }

    memset(&rebuilt, 0, sizeof(rebuilt));
    input_points = (TMPoint *) calloc(mesh->point_count + 1, sizeof(*input_points));
    if (input_points == NULL) {
        return TM_ERR_ALLOC;
    }

    memcpy(input_points, mesh->points, mesh->point_count * sizeof(*input_points));
    point_index = (int) mesh->point_count;
    input_points[point_index].xy[0] = split_point[0];
    input_points[point_index].xy[1] = split_point[1];
    input_points[point_index].original_index = point_index;
    input_points[point_index].kind = TM_VERTEX_SEGMENT_SPLIT;
    input_points[point_index].incident_triangle = -1;
    input_points[point_index].protection_apex = -1;
    input_points[point_index].protection_side = TM_PROTECTION_SIDE_NONE;
    input_points[point_index].protection_level = 0;

    status = tm_build_mesh(input_points, mesh->point_count + 1, &rebuilt);
    if (status == TM_OK) {
        size_t i;

        tm_verbose_log(
            options,
            "recover segment %d-%d: rebuilt unconstrained triangulation with %zu points and %zu triangles",
            a,
            b,
            rebuilt.point_count,
            rebuilt.triangle_count
        );
        for (i = 0; i < rebuilt.point_count; ++i) {
            rebuilt.points[i].kind = input_points[i].kind;
            rebuilt.points[i].original_index = input_points[i].original_index;
            rebuilt.points[i].protection_apex = input_points[i].protection_apex;
            rebuilt.points[i].protection_side = input_points[i].protection_side;
            rebuilt.points[i].protection_level = input_points[i].protection_level;
        }

        status = tm_copy_live_segments_with_split(mesh, &rebuilt, a, b, point_index);
        if (status != TM_OK) {
            tm_verbose_log(options, "recover segment %d-%d: copying live segments to rebuilt mesh failed", a, b);
        }
    } else {
        tm_verbose_log(options, "recover segment %d-%d: rebuilding unconstrained triangulation failed", a, b);
    }
    if (status == TM_OK) {
        tm_verbose_log(options, "recover segment %d-%d: recovering constraints on rebuilt mesh", a, b);
        status = tm_recover_all_live_segments(&rebuilt, options);
        if (status != TM_OK) {
            tm_verbose_log(options, "recover segment %d-%d: constraint recovery on rebuilt mesh failed", a, b);
        }
    }

    free(input_points);
    if (status != TM_OK) {
        tm_verbose_log(options, "recover segment %d-%d: rebuilt recovery path failed", a, b);
        tm_free_mesh(&rebuilt);
        return status;
    }

    tm_free_mesh(mesh);
    *mesh = rebuilt;
    return status;
}

static TMStatus tm_recover_segment(TMMesh *mesh, int a, int b, const TMBuildOptions *options)
{
    size_t iteration_limit;
    size_t iteration;

    if (tm_mesh_has_edge(mesh, a, b, NULL, NULL)) {
        return tm_mark_constraint_edge(mesh, a, b);
    }

    iteration_limit = mesh->triangle_count * mesh->triangle_count + mesh->point_count * 8 + 16;
    for (iteration = 0; iteration < iteration_limit; ++iteration) {
        int triangle_index;
        int edge;
        TMStatus status;

        status = tm_find_crossing_edge(mesh, a, b, &triangle_index, &edge);
        if (status != TM_OK) {
            break;
        }

        tm_verbose_log(options, "recover segment %d-%d: flip crossing edge in triangle %d edge %d", a, b, triangle_index, edge);
        status = tm_flip_edge(mesh, triangle_index, edge);
        if (status != TM_OK) {
            return status;
        }

        if (tm_mesh_has_edge(mesh, a, b, NULL, NULL)) {
            return tm_mark_constraint_edge(mesh, a, b);
        }
    }

    if (tm_mesh_has_edge(mesh, a, b, NULL, NULL)) {
        return tm_mark_constraint_edge(mesh, a, b);
    }

    return tm_split_recovery_segment(mesh, a, b, options);
}

static TMStatus tm_restore_constrained_delaunay(TMMesh *mesh, const TMBuildOptions *options)
{
    int changed = 1;

    while (changed) {
        size_t tri_index;

        changed = 0;
        for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
            int edge;

            for (edge = 0; edge < 3; ++edge) {
                int neighbor = mesh->triangles[tri_index].nbr[edge];

                if (neighbor < 0 || (size_t) neighbor <= tri_index || mesh->triangles[tri_index].constrained[edge]) {
                    continue;
                }

                if (!tm_edge_is_locally_delaunay(mesh, (int) tri_index, edge)) {
                    TMStatus status;

                    tm_verbose_log(options, "restore CDT: flip edge in triangle %zu edge %d", tri_index, edge);
                    status = tm_flip_edge(mesh, (int) tri_index, edge);
                    if (status != TM_OK) {
                        return status;
                    }
                    changed = 1;
                    break;
                }
            }

            if (changed) {
                break;
            }
        }
    }

    return TM_OK;
}

static TMStatus tm_update_after_local_edit(TMMesh *mesh, const TMBuildOptions *options)
{
    TMStatus status = tm_mark_all_constraints(mesh);

    if (status != TM_OK) {
        return status;
    }

    return tm_restore_constrained_delaunay(mesh, options);
}

static TMStatus tm_split_segment_with_point(
    TMMesh *mesh,
    size_t segment_index,
    const double point[2],
    TMVertexKind kind,
    const TMBuildOptions *options,
    int *out_point_index
)
{
    const TMSegment *segment = &mesh->segments[segment_index];
    int triangle_index;
    int edge;
    int point_index = -1;
    TMStatus status;

    if (!tm_mesh_has_edge(mesh, segment->v[0], segment->v[1], &triangle_index, &edge)) {
        return TM_ERR_INVALID_MESH;
    }

    status = tm_insert_point_on_edge(mesh, triangle_index, edge, point, kind, &point_index);
    if (status != TM_OK) {
        return status;
    }

    status = tm_split_segment_registry(mesh, segment_index, point_index);
    if (status != TM_OK) {
        return status;
    }

    status = tm_update_after_local_edit(mesh, options);
    if (status != TM_OK) {
        return status;
    }

    if (out_point_index != NULL) {
        *out_point_index = point_index;
    }
    return TM_OK;
}

static TMStatus tm_split_segment_midpoint(
    TMMesh *mesh,
    size_t segment_index,
    const TMBuildOptions *options,
    int *out_point_index
)
{
    double midpoint[2];
    const TMSegment *segment = &mesh->segments[segment_index];

    midpoint[0] = 0.5 * (mesh->points[segment->v[0]].xy[0] + mesh->points[segment->v[1]].xy[0]);
    midpoint[1] = 0.5 * (mesh->points[segment->v[0]].xy[1] + mesh->points[segment->v[1]].xy[1]);
    return tm_split_segment_with_point(mesh, segment_index, midpoint, TM_VERTEX_SEGMENT_SPLIT, options, out_point_index);
}

static TMStatus tm_mark_live_segment_protected_with_apex(TMMesh *mesh, int a, int b, int apex)
{
    size_t segment_index;

    if (!tm_find_live_segment_index(mesh, a, b, &segment_index)) {
        return TM_ERR_INVALID_MESH;
    }

    mesh->segments[segment_index].is_protected = 1;
    mesh->segments[segment_index].protected_apex = apex;
    return TM_OK;
}

static TMStatus tm_mark_live_segment_protected(TMMesh *mesh, int apex, int other)
{
    return tm_mark_live_segment_protected_with_apex(mesh, apex, other, apex);
}

static void tm_mark_shell_point(TMMesh *mesh, int point_index, int apex, TMProtectionSide side, unsigned int level)
{
    mesh->points[point_index].protection_apex = apex;
    mesh->points[point_index].protection_side = (unsigned char) side;
    mesh->points[point_index].protection_level = level;
}

static size_t tm_local_count_protected_corners(const TMMesh *mesh)
{
    return tm_count_protected_corners(mesh);
}

static size_t tm_count_exempt_triangles(const TMMesh *mesh)
{
    size_t tri_index;
    size_t count = 0;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        if (tm_local_triangle_is_exempt_from_quality(mesh, (int) tri_index)) {
            count += 1;
        }
    }

    return count;
}

static TMStatus tm_apply_simple_acute_corner_protection(
    TMMesh *mesh,
    size_t incident_segments[2],
    int apex,
    size_t levels,
    const TMBuildOptions *options
)
{
    size_t level;
    int left_point = -1;
    int right_point = -1;
    TMStatus status;

    for (level = 0; level < levels; ++level) {
        status = tm_split_segment_midpoint(mesh, incident_segments[0], options, &left_point);
        if (status != TM_OK) {
            return status;
        }
        if (!tm_find_live_segment_index(mesh, apex, left_point, &incident_segments[0])) {
            return TM_ERR_INVALID_MESH;
        }

        status = tm_split_segment_midpoint(mesh, incident_segments[1], options, &right_point);
        if (status != TM_OK) {
            return status;
        }
        if (!tm_find_live_segment_index(mesh, apex, right_point, &incident_segments[1])) {
            return TM_ERR_INVALID_MESH;
        }
    }

    status = tm_mark_live_segment_protected(mesh, apex, left_point);
    if (status == TM_OK) {
        status = tm_mark_live_segment_protected(mesh, apex, right_point);
    }
    return status;
}

static TMStatus tm_build_shell_chain(
    TMMesh *mesh,
    int apex,
    int outer,
    size_t segment_index,
    TMProtectionSide side,
    double base_length,
    size_t levels,
    const TMBuildOptions *options
)
{
    double segment_length = sqrt(tm_segment_length_squared(mesh, apex, outer));
    double direction[2];
    int inner = apex;
    size_t level;

    if (segment_length <= 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    direction[0] = (mesh->points[outer].xy[0] - mesh->points[apex].xy[0]) / segment_length;
    direction[1] = (mesh->points[outer].xy[1] - mesh->points[apex].xy[1]) / segment_length;

    for (level = 1; level <= levels; ++level) {
        double distance = base_length * (double) (1ULL << (level - 1));
        double point[2];
        int point_index;
        TMStatus status;

        if (!(distance < segment_length * (1.0 - 1e-12))) {
            return TM_ERR_INVALID_MESH;
        }

        point[0] = mesh->points[apex].xy[0] + direction[0] * distance;
        point[1] = mesh->points[apex].xy[1] + direction[1] * distance;

        status = tm_split_segment_with_point(mesh, segment_index, point, TM_VERTEX_SEGMENT_SPLIT, options, &point_index);
        if (status != TM_OK) {
            return status;
        }

        tm_mark_shell_point(mesh, point_index, apex, side, (unsigned int) level);
        status = tm_mark_live_segment_protected_with_apex(mesh, inner, point_index, apex);
        if (status != TM_OK) {
            return status;
        }

        inner = point_index;
        if (level < levels && !tm_find_live_segment_index(mesh, inner, outer, &segment_index)) {
            return TM_ERR_INVALID_MESH;
        }
    }

    return TM_OK;
}

static TMStatus tm_apply_shell_acute_corner_protection(
    TMMesh *mesh,
    size_t incident_segments[2],
    int apex,
    int left,
    int right,
    size_t levels,
    const TMBuildOptions *options
)
{
    double left_length = sqrt(tm_segment_length_squared(mesh, apex, left));
    double right_length = sqrt(tm_segment_length_squared(mesh, apex, right));
    double base_length = (left_length < right_length) ? left_length : right_length;
    TMStatus status;

    base_length /= (double) (1ULL << levels);
    if (base_length <= 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    status = tm_build_shell_chain(
        mesh,
        apex,
        left,
        incident_segments[0],
        TM_PROTECTION_SIDE_LEFT,
        base_length,
        levels,
        options
    );
    if (status != TM_OK) {
        return status;
    }

    return tm_build_shell_chain(
        mesh,
        apex,
        right,
        incident_segments[1],
        TM_PROTECTION_SIDE_RIGHT,
        base_length,
        levels,
        options
    );
}

static TMStatus tm_apply_acute_corner_protection(TMMesh *mesh, const TMBuildOptions *options)
{
    double protect_angle_deg = tm_protection_angle_deg(options);
    size_t protection_level_limit = tm_protection_level_limit(options);
    int apex;

    if (!tm_protect_acute_corners(options)) {
        return TM_OK;
    }

    for (apex = 0; (size_t) apex < mesh->point_count; ++apex) {
        size_t incident_segments[2];
        int left;
        int right;
        double angle_deg;
        size_t levels = 1;
        TMStatus status;

        if (mesh->points[apex].kind != TM_VERTEX_INPUT) {
            continue;
        }
        if (!tm_find_live_incident_segments(mesh, apex, incident_segments)) {
            continue;
        }
        if (mesh->segments[incident_segments[0]].is_protected || mesh->segments[incident_segments[1]].is_protected) {
            continue;
        }

        left = tm_segment_other_endpoint(&mesh->segments[incident_segments[0]], apex);
        right = tm_segment_other_endpoint(&mesh->segments[incident_segments[1]], apex);
        angle_deg = tm_corner_angle_deg(mesh, apex, left, right);
        if (!(angle_deg < protect_angle_deg * (1.0 - 1e-12))) {
            continue;
        }

        while (levels < protection_level_limit && angle_deg * (double) (1ULL << levels) < protect_angle_deg) {
            levels += 1;
        }

        tm_verbose_log(
            options,
            "acute protection: apex %d angle %.6f deg, splitting %zu level(s) in %s mode",
            apex,
            angle_deg,
            levels,
            tm_acute_mode_name(options)
        );

        if (tm_acute_mode(options) == TM_ACUTE_MODE_SIMPLE) {
            status = tm_apply_simple_acute_corner_protection(mesh, incident_segments, apex, levels, options);
        } else {
            status = tm_apply_shell_acute_corner_protection(mesh, incident_segments, apex, left, right, levels, options);
        }
        if (status != TM_OK) {
            return status;
        }
    }

    if (tm_local_count_protected_corners(mesh) != 0) {
        tm_verbose_log(
            options,
            "acute protection ready: %zu protected corners, %zu exempt shield triangles",
            tm_local_count_protected_corners(mesh),
            tm_count_exempt_triangles(mesh)
        );
    }

    return TM_OK;
}

static TMStatus tm_insert_triangle_split_point(
    TMMesh *mesh,
    int start_triangle,
    const double point[2],
    const TMBuildOptions *options,
    int *out_point_index
)
{
    TMLocation location;
    int on_edge = -1;
    TMStatus status;

    status = tm_locate_point(mesh, point, start_triangle, &location);
    if (status != TM_OK) {
        return status;
    }

    if (!tm_point_inside_or_on_triangle(mesh, location.triangle, point, &on_edge)) {
        return TM_ERR_INVALID_MESH;
    }

    if (on_edge >= 0) {
        int a;
        int b;
        size_t segment_index;

        if (mesh->triangles[location.triangle].constrained[on_edge]) {
            tm_triangle_edge_vertices(&mesh->triangles[location.triangle], on_edge, &a, &b);
            if (!tm_find_live_segment_index(mesh, a, b, &segment_index)) {
                return TM_ERR_INVALID_MESH;
            }
            return tm_split_segment_midpoint(mesh, segment_index, options, out_point_index);
        }

        status = tm_insert_point_on_edge(
            mesh,
            location.triangle,
            on_edge,
            point,
            TM_VERTEX_TRIANGLE_SPLIT,
            out_point_index
        );
    } else {
        status = tm_insert_point_in_triangle(
            mesh,
            location.triangle,
            point,
            TM_VERTEX_TRIANGLE_SPLIT,
            out_point_index
        );
    }

    if (status != TM_OK) {
        return status;
    }

    return tm_update_after_local_edit(mesh, options);
}

static TMStatus tm_refine_quality_mesh(TMMesh *mesh, const TMBuildOptions *options)
{
    double min_angle_deg = tm_quality_min_angle_deg(options);
    double max_area = tm_area_limit(options);
    double beta = tm_quality_beta(min_angle_deg);
    size_t step_limit = tm_refinement_step_limit(options, mesh);
    size_t steps = 0;
    size_t segment_splits = 0;
    size_t triangle_splits = 0;

    if (max_area > 0.0) {
        tm_verbose_log(
            options,
            "quality refinement: target min angle %.2f deg (beta=%.6f), max area=%.6f, step limit=%zu, mode=%s",
            min_angle_deg,
            beta,
            max_area,
            step_limit,
            tm_use_offcenters(options) ? "off-centers" : "circumcenters"
        );
    } else {
        tm_verbose_log(
            options,
            "quality refinement: target min angle %.2f deg (beta=%.6f), max area=off, step limit=%zu, mode=%s",
            min_angle_deg,
            beta,
            step_limit,
            tm_use_offcenters(options) ? "off-centers" : "circumcenters"
        );
    }

    while (steps < step_limit) {
        size_t segment_index;
        int encroaching_point = -1;
        size_t triangle_index;
        double ratio;
        double area;
        int quality_violation = 0;
        TMStatus status;

        if (tm_find_encroached_segment(mesh, &segment_index, &encroaching_point)) {
            const TMSegment *segment = &mesh->segments[segment_index];

            tm_verbose_log(
                options,
                "refine step %zu: split encroached segment %d-%d due to point %d",
                steps + 1,
                segment->v[0],
                segment->v[1],
                encroaching_point
            );
            status = tm_split_segment_midpoint(mesh, segment_index, options, NULL);
            if (status != TM_OK) {
                return status;
            }
            segment_splits += 1;
            steps += 1;
            continue;
        }

        if (!tm_find_bad_triangle(
                mesh,
                beta,
                max_area,
                tm_use_offcenters(options),
                &triangle_index,
                &ratio,
                &area,
                &quality_violation
            )) {
            tm_verbose_log(
                options,
                "quality refinement complete: %zu steps, %zu segment splits, %zu triangle splits",
                steps,
                segment_splits,
                triangle_splits
            );
            return TM_OK;
        }

        {
            double circumcenter[2];
            double offcenter[2];
            double candidate[2];
            const char *candidate_name = "circumcenter";
            size_t circumcenter_segment_index = 0;
            size_t offcenter_segment_index = 0;
            size_t candidate_segment_index = 0;
            size_t split_segment_index = 0;
            int circumcenter_encroaches = 0;
            int candidate_encroaches = 0;
            int offcenter_available = 0;
            int offcenter_encroaches = 0;
            int split_segment = 0;
            int use_offcenter = quality_violation && tm_use_offcenters(options);

            status = tm_triangle_circumcenter(mesh, (int) triangle_index, circumcenter);
            if (status != TM_OK) {
                return status;
            }

            circumcenter_encroaches = tm_find_encroached_segment_for_point(mesh, circumcenter, &circumcenter_segment_index);
            candidate[0] = circumcenter[0];
            candidate[1] = circumcenter[1];
            candidate_encroaches = circumcenter_encroaches;
            candidate_segment_index = circumcenter_segment_index;

            if (use_offcenter) {
                offcenter_available = tm_triangle_offcenter_candidate(
                    mesh,
                    (int) triangle_index,
                    min_angle_deg,
                    circumcenter,
                    offcenter
                );
                if (offcenter_available) {
                    offcenter_encroaches = tm_find_encroached_segment_for_point(mesh, offcenter, &offcenter_segment_index);
                    if (!offcenter_encroaches) {
                        candidate[0] = offcenter[0];
                        candidate[1] = offcenter[1];
                        candidate_name = "off-center";
                        candidate_encroaches = 0;
                    } else if (circumcenter_encroaches) {
                        double offcenter_length_sq = tm_segment_length_squared(
                            mesh,
                            mesh->segments[offcenter_segment_index].v[0],
                            mesh->segments[offcenter_segment_index].v[1]
                        );
                        double circumcenter_length_sq = tm_segment_length_squared(
                            mesh,
                            mesh->segments[circumcenter_segment_index].v[0],
                            mesh->segments[circumcenter_segment_index].v[1]
                        );

                        split_segment = 1;
                        split_segment_index = (offcenter_length_sq < circumcenter_length_sq) ?
                            offcenter_segment_index :
                            circumcenter_segment_index;
                    }
                }
            }

            if (!split_segment && candidate_encroaches) {
                split_segment = 1;
                split_segment_index = candidate_segment_index;
            }

            if (split_segment) {
                const TMSegment *segment = &mesh->segments[split_segment_index];

                tm_verbose_log(
                    options,
                    quality_violation ?
                        "refine step %zu: bad triangle %zu (R/s=%.6f, area=%.6f) encroaches segment %d-%d, split segment first" :
                        "refine step %zu: oversized triangle %zu (area=%.6f, R/s=%.6f) encroaches segment %d-%d, split segment first",
                    steps + 1,
                    triangle_index,
                    quality_violation ? ratio : area,
                    quality_violation ? area : ratio,
                    segment->v[0],
                    segment->v[1]
                );
                status = tm_split_segment_midpoint(mesh, split_segment_index, options, NULL);
                if (status != TM_OK) {
                    return status;
                }
                segment_splits += 1;
                steps += 1;
                continue;
            }

            tm_verbose_log(
                options,
                quality_violation ?
                    "refine step %zu: split bad triangle %zu at %s (R/s=%.6f, area=%.6f)" :
                    "refine step %zu: split oversized triangle %zu at %s (area=%.6f, R/s=%.6f)",
                steps + 1,
                triangle_index,
                candidate_name,
                quality_violation ? ratio : area,
                quality_violation ? area : ratio
            );
            status = tm_insert_triangle_split_point(mesh, (int) triangle_index, candidate, options, NULL);
            if (status != TM_OK) {
                return status;
            }
            triangle_splits += 1;
            steps += 1;
        }
    }

    return TM_ERR_INTERNAL;
}

static int tm_point_in_domain(const TMPSLG *pslg, const double point[2])
{
    size_t segment_index;
    int winding = 0;

    for (segment_index = 0; segment_index < pslg->segment_count; ++segment_index) {
        const double *a = pslg->points[pslg->segments[segment_index].v[0]].xy;
        const double *b = pslg->points[pslg->segments[segment_index].v[1]].xy;
        double ay = a[1];
        double by = b[1];

        if ((ay > point[1]) == (by > point[1])) {
            continue;
        }

        {
            double x = a[0] + (point[1] - ay) * (b[0] - a[0]) / (by - ay);
            if (x > point[0]) {
                winding = !winding;
            }
        }
    }

    return winding;
}

static TMStatus tm_filter_domain_triangles(TMMesh *mesh, const TMPSLG *pslg, const TMBuildOptions *options)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        double centroid[2];
        const TMTriangle *triangle = &mesh->triangles[tri_index];

        centroid[0] = (mesh->points[triangle->v[0]].xy[0] +
                       mesh->points[triangle->v[1]].xy[0] +
                       mesh->points[triangle->v[2]].xy[0]) / 3.0;
        centroid[1] = (mesh->points[triangle->v[0]].xy[1] +
                       mesh->points[triangle->v[1]].xy[1] +
                       mesh->points[triangle->v[2]].xy[1]) / 3.0;

        if (!tm_point_in_domain(pslg, centroid)) {
            mesh->triangles[tri_index].dead = 1;
        }
    }

    tm_verbose_log(options, "domain filter: compacting kept triangles");

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

    if (mesh->triangle_count == 0) {
        return TM_ERR_INVALID_MESH;
    }

    if (tm_rebuild_topology(mesh) != TM_OK) {
        return TM_ERR_INVALID_MESH;
    }

    return tm_mark_all_constraints(mesh);
}

static TMStatus tm_assign_coverage_regions(
    const TMMesh *mesh,
    const TMRegion *regions,
    size_t region_count,
    int **out_triangle_markers
)
{
    int *triangle_markers = NULL;
    int *stack = NULL;
    size_t region_index;

    if (mesh == NULL || out_triangle_markers == NULL || region_count == 0 || regions == NULL) {
        return TM_ERR_INTERNAL;
    }

    *out_triangle_markers = NULL;
    triangle_markers = (int *) malloc(mesh->triangle_count * sizeof(*triangle_markers));
    stack = (int *) malloc(mesh->triangle_count * sizeof(*stack));
    if (triangle_markers == NULL || stack == NULL) {
        free(triangle_markers);
        free(stack);
        return TM_ERR_ALLOC;
    }

    for (region_index = 0; region_index < mesh->triangle_count; ++region_index) {
        triangle_markers[region_index] = INT_MIN;
    }

    for (region_index = 0; region_index < region_count; ++region_index) {
        TMLocation location;
        int seed_triangle;
        size_t stack_size = 0;
        TMStatus status = tm_locate_point(mesh, regions[region_index].xy, -1, &location);

        if (status != TM_OK || location.triangle < 0 || (size_t) location.triangle >= mesh->triangle_count) {
            tm_set_pslg_error_detail(
                "coverage region seed %zu at (%.12g, %.12g) could not be located inside the triangulation",
                region_index,
                regions[region_index].xy[0],
                regions[region_index].xy[1]
            );
            free(triangle_markers);
            free(stack);
            return TM_ERR_INVALID_PSLG;
        }

        seed_triangle = location.triangle;
        if (triangle_markers[seed_triangle] != INT_MIN &&
            triangle_markers[seed_triangle] != regions[region_index].marker) {
            tm_set_pslg_error_detail(
                "coverage region seed %zu at (%.12g, %.12g) conflicts with marker %d already assigned to triangle %d",
                region_index,
                regions[region_index].xy[0],
                regions[region_index].xy[1],
                triangle_markers[seed_triangle],
                seed_triangle
            );
            free(triangle_markers);
            free(stack);
            return TM_ERR_INVALID_PSLG;
        }

        triangle_markers[seed_triangle] = regions[region_index].marker;
        stack[stack_size++] = seed_triangle;
        while (stack_size > 0) {
            int triangle_index = stack[--stack_size];
            const TMTriangle *triangle = &mesh->triangles[triangle_index];
            int edge;

            if (triangle_markers[triangle_index] != INT_MIN &&
                triangle_markers[triangle_index] != regions[region_index].marker) {
                tm_set_pslg_error_detail(
                    "coverage region marker %d leaks into triangle %d already claimed by marker %d",
                    regions[region_index].marker,
                    triangle_index,
                    triangle_markers[triangle_index]
                );
                free(triangle_markers);
                free(stack);
                return TM_ERR_INVALID_PSLG;
            }

            for (edge = 0; edge < 3; ++edge) {
                int neighbor = triangle->nbr[edge];

                if (triangle->constrained[edge] || neighbor < 0) {
                    continue;
                }
                if (triangle_markers[neighbor] == INT_MIN) {
                    triangle_markers[neighbor] = regions[region_index].marker;
                    stack[stack_size++] = neighbor;
                    continue;
                }
                if (triangle_markers[neighbor] != regions[region_index].marker) {
                    tm_set_pslg_error_detail(
                        "coverage region marker %d leaks into triangle %d already claimed by marker %d",
                        regions[region_index].marker,
                        neighbor,
                        triangle_markers[neighbor]
                    );
                    free(triangle_markers);
                    free(stack);
                    return TM_ERR_INVALID_PSLG;
                }
            }
        }
    }

    *out_triangle_markers = triangle_markers;
    free(stack);
    return TM_OK;
}

static TMStatus tm_compact_marked_triangles(TMMesh *mesh, int **triangle_markers)
{
    size_t read_index;
    size_t write_index = 0;
    int *markers = NULL;

    if (mesh == NULL || triangle_markers == NULL || *triangle_markers == NULL) {
        return TM_ERR_INTERNAL;
    }

    markers = *triangle_markers;
    for (read_index = 0; read_index < mesh->triangle_count; ++read_index) {
        if (markers[read_index] == INT_MIN) {
            continue;
        }
        if (write_index != read_index) {
            mesh->triangles[write_index] = mesh->triangles[read_index];
            markers[write_index] = markers[read_index];
        }
        mesh->triangles[write_index].dead = 0;
        write_index += 1;
    }

    mesh->triangle_count = write_index;
    if (mesh->triangle_count == 0) {
        tm_set_pslg_error_detail("coverage graph produced no marked triangles");
        return TM_ERR_INVALID_MESH;
    }

    if (tm_rebuild_topology(mesh) != TM_OK) {
        return TM_ERR_INVALID_MESH;
    }

    return tm_mark_all_constraints(mesh);
}

TMStatus tm_build_pslg_mesh(const TMPSLG *pslg, const TMBuildOptions *options, TMMesh *out_mesh)
{
    TMPSLG working_pslg;
    size_t segment_index;
    TMStatus status;

    if (pslg == NULL || out_mesh == NULL) {
        return TM_ERR_INTERNAL;
    }

    memset(&working_pslg, 0, sizeof(working_pslg));
    tm_verbose_log(options, "read PSLG: %zu vertices, %zu segments, %zu hole markers", pslg->point_count, pslg->segment_count, pslg->hole_count);

    status = tm_clone_pslg(pslg, &working_pslg);
    if (status != TM_OK) {
        return status;
    }

    status = tm_pre_split_long_segments(&working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        return status;
    }
    status = tm_pre_split_near_segment_free_points(&working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        return status;
    }
    status = tm_seed_interior_pslg_points(&working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        return status;
    }
    tm_verbose_log(
        options,
        "PSLG after blocker pre-split: %zu vertices, %zu segments",
        working_pslg.point_count,
        working_pslg.segment_count
    );

    status = tm_build_mesh(working_pslg.points, working_pslg.point_count, out_mesh);
    if (status != TM_OK) {
        tm_verbose_log(options, "failed to build initial triangulation from pre-split PSLG");
        tm_free_pslg(&working_pslg);
        return status;
    }

    status = tm_copy_segments_to_mesh(&working_pslg, out_mesh);
    if (status != TM_OK) {
        tm_verbose_log(options, "failed to attach pre-split segments to mesh");
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    tm_verbose_log(options, "initial triangulation: %zu triangles", out_mesh->triangle_count);

    for (segment_index = 0; segment_index < out_mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &out_mesh->segments[segment_index];

        tm_verbose_log(
            options,
            "recovering segment %zu/%zu: %d-%d",
            segment_index + 1,
            out_mesh->segment_count,
            segment->v[0],
            segment->v[1]
        );

        status = tm_recover_segment(out_mesh, segment->v[0], segment->v[1], options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }
    }

    tm_verbose_log(options, "restoring constrained Delaunay condition");
    status = tm_restore_constrained_delaunay(out_mesh, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    tm_verbose_log(options, "classifying triangles against PSLG domain");
    status = tm_filter_domain_triangles(out_mesh, &working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    if (tm_refinement_enabled(options)) {
        status = tm_apply_acute_corner_protection(out_mesh, options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }

        status = tm_refine_quality_mesh(out_mesh, options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }
    }

    tm_verbose_log(
        options,
        "final mesh: %zu points, %zu triangles, %zu live subsegments",
        out_mesh->point_count,
        out_mesh->triangle_count,
        tm_count_live_segments(out_mesh)
    );
    if (tm_local_count_protected_corners(out_mesh) != 0) {
        tm_verbose_log(
            options,
            "protected corners: %zu, exempt triangles: %zu",
            tm_local_count_protected_corners(out_mesh),
            tm_count_exempt_triangles(out_mesh)
        );
    }
    tm_free_pslg(&working_pslg);
    return TM_OK;
}

TMStatus tm_build_coverage_mesh(
    const TMPSLG *pslg,
    const TMRegion *regions,
    size_t region_count,
    const TMBuildOptions *options,
    TMMesh *out_mesh,
    int **out_triangle_markers
)
{
    TMPSLG working_pslg;
    size_t segment_index;
    TMStatus status;

    if (pslg == NULL || out_mesh == NULL || out_triangle_markers == NULL) {
        return TM_ERR_INTERNAL;
    }
    if (regions == NULL || region_count == 0) {
        tm_set_pslg_error_detail("coverage meshing requires at least one region seed");
        return TM_ERR_INVALID_PSLG;
    }

    *out_triangle_markers = NULL;
    memset(&working_pslg, 0, sizeof(working_pslg));
    tm_verbose_log(
        options,
        "read coverage graph: %zu vertices, %zu segments, %zu regions",
        pslg->point_count,
        pslg->segment_count,
        region_count
    );

    status = tm_validate_segment_graph(pslg->points, pslg->point_count, pslg->segments, pslg->segment_count);
    if (status != TM_OK) {
        return status;
    }

    status = tm_clone_pslg(pslg, &working_pslg);
    if (status != TM_OK) {
        return status;
    }

    status = tm_pre_split_long_segments(&working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        return status;
    }
    status = tm_pre_split_near_segment_free_points(&working_pslg, options);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        return status;
    }
    tm_verbose_log(
        options,
        "coverage graph after blocker pre-split: %zu vertices, %zu segments",
        working_pslg.point_count,
        working_pslg.segment_count
    );

    status = tm_build_mesh(working_pslg.points, working_pslg.point_count, out_mesh);
    if (status != TM_OK) {
        tm_verbose_log(options, "failed to build initial triangulation from coverage graph");
        tm_free_pslg(&working_pslg);
        return status;
    }

    status = tm_copy_segments_to_mesh(&working_pslg, out_mesh);
    if (status != TM_OK) {
        tm_verbose_log(options, "failed to attach coverage segments to mesh");
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    tm_verbose_log(options, "initial coverage triangulation: %zu triangles", out_mesh->triangle_count);

    for (segment_index = 0; segment_index < out_mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &out_mesh->segments[segment_index];

        tm_verbose_log(
            options,
            "recovering coverage segment %zu/%zu: %d-%d",
            segment_index + 1,
            out_mesh->segment_count,
            segment->v[0],
            segment->v[1]
        );

        status = tm_recover_segment(out_mesh, segment->v[0], segment->v[1], options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }
    }

    tm_verbose_log(options, "marking constrained edges for coverage mesh");
    status = tm_mark_all_constraints(out_mesh);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    if (tm_refinement_enabled(options)) {
        status = tm_apply_acute_corner_protection(out_mesh, options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }

        status = tm_refine_quality_mesh(out_mesh, options);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }
    }

    status = tm_assign_coverage_regions(out_mesh, regions, region_count, out_triangle_markers);
    if (status != TM_OK) {
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    status = tm_compact_marked_triangles(out_mesh, out_triangle_markers);
    if (status != TM_OK) {
        free(*out_triangle_markers);
        *out_triangle_markers = NULL;
        tm_free_pslg(&working_pslg);
        tm_free_mesh(out_mesh);
        return status;
    }

    tm_verbose_log(
        options,
        "final coverage mesh: %zu points, %zu triangles, %zu live subsegments",
        out_mesh->point_count,
        out_mesh->triangle_count,
        tm_count_live_segments(out_mesh)
    );
    tm_free_pslg(&working_pslg);
    return TM_OK;
}
