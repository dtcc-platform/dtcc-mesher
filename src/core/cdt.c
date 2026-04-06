#include "cdt.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>

#include "predicates.h"
#include "validate.h"

static const double tm_pi = 3.14159265358979323846;
static const double tm_default_min_angle_deg = 20.0;
static const double tm_default_protection_angle_deg = 60.0;
static const size_t tm_default_max_protection_levels = 24;

typedef struct {
    size_t segment_index;
    int point_index;
} TMEncroachmentCandidate;

typedef struct {
    size_t triangle_index;
} TMBadTriangleCandidate;

static int tm_find_live_segment_index(const TMMesh *mesh, int a, int b, size_t *out_segment_index);
static TMStatus tm_append_live_segment_unique(
    TMMesh *mesh,
    int a,
    int b,
    int original_index,
    const TMSegment *source_flags
);
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
static TMStatus tm_restore_constrained_delaunay(TMMesh *mesh, const TMBuildOptions *options);
static int tm_evaluate_triangle_refinement_need(
    const TMMesh *mesh,
    size_t triangle_index,
    double beta,
    double max_area,
    double *out_ratio,
    double *out_area,
    int *out_quality_violation
);

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

static double tm_elapsed_seconds(clock_t start, clock_t end)
{
    return (double) (end - start) / (double) CLOCKS_PER_SEC;
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

    return tm_default_protection_angle_deg;
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

    return tm_default_max_protection_levels;
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
    double scale;
    double tolerance;
    double t;
    double projection[2];
    double distance_sq;

    if (len_sq == 0.0) {
        return 0;
    }

    scale = fabs(a[0]) + fabs(a[1]) + fabs(b[0]) + fabs(b[1]) + fabs(point[0]) + fabs(point[1]);
    tolerance = 64.0 * DBL_EPSILON * ((scale > 1.0) ? scale : 1.0);

    t = ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / len_sq;
    if (t < -tolerance || t > 1.0 + tolerance) {
        return 0;
    }

    projection[0] = a[0] + t * dx;
    projection[1] = a[1] + t * dy;
    distance_sq = (point[0] - projection[0]) * (point[0] - projection[0]) +
                  (point[1] - projection[1]) * (point[1] - projection[1]);
    return distance_sq <= tolerance * tolerance;
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

static double tm_segment_split_tolerance_xy(const double a[2], const double b[2], const double point[2])
{
    double scale = fabs(a[0]) + fabs(a[1]) + fabs(b[0]) + fabs(b[1]) + fabs(point[0]) + fabs(point[1]);

    return 64.0 * DBL_EPSILON * ((scale > 1.0) ? scale : 1.0);
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
    unsigned char *flags = NULL;
    size_t tri_index;
    size_t segment_index;

    if (mesh->triangle_count != 0) {
        flags = (unsigned char *) calloc(mesh->triangle_count * 3u, sizeof(*flags));
        if (flags == NULL) {
            return TM_ERR_ALLOC;
        }
    }

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        int found = 0;

        if (!mesh->segments[segment_index].live) {
            continue;
        }
        for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
            int edge = tm_find_edge_in_triangle(
                &mesh->triangles[tri_index],
                mesh->segments[segment_index].v[0],
                mesh->segments[segment_index].v[1]
            );

            if (edge >= 0) {
                flags[tri_index * 3u + (size_t) edge] = 1;
                found = 1;
            }
        }
        if (!found) {
            free(flags);
            tm_set_pslg_error_detail(
                "live segment %d-%d is not present as a mesh edge after constraint recovery",
                mesh->segments[segment_index].v[0],
                mesh->segments[segment_index].v[1]
            );
            return TM_ERR_INVALID_MESH;
        }
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            mesh->triangles[tri_index].constrained[edge] = flags[tri_index * 3u + (size_t) edge];
        }
    }

    free(flags);

    return TM_OK;
}

static const TMSegment *tm_find_snapshot_segment(
    const TMSegment *segments,
    size_t segment_count,
    int a,
    int b
)
{
    size_t segment_index;
    int lo = (a < b) ? a : b;
    int hi = (a < b) ? b : a;

    for (segment_index = 0; segment_index < segment_count; ++segment_index) {
        int sa;
        int sb;
        int slo;
        int shi;

        if (!segments[segment_index].live) {
            continue;
        }

        sa = segments[segment_index].v[0];
        sb = segments[segment_index].v[1];
        slo = (sa < sb) ? sa : sb;
        shi = (sa < sb) ? sb : sa;
        if (slo == lo && shi == hi) {
            return &segments[segment_index];
        }
    }

    return NULL;
}

static TMStatus tm_sync_live_segments_from_constraints(TMMesh *mesh)
{
    TMSegment *snapshot = NULL;
    size_t snapshot_count = 0;
    size_t tri_index;

    if (mesh == NULL) {
        return TM_ERR_INTERNAL;
    }

    if (mesh->segment_count != 0) {
        snapshot = (TMSegment *) malloc(mesh->segment_count * sizeof(*snapshot));
        if (snapshot == NULL) {
            return TM_ERR_ALLOC;
        }
        memcpy(snapshot, mesh->segments, mesh->segment_count * sizeof(*snapshot));
        snapshot_count = mesh->segment_count;
    }

    mesh->segment_count = 0;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            int a;
            int b;
            const TMSegment *source_segment;
            int original_index = -1;
            TMStatus status;

            if (!mesh->triangles[tri_index].constrained[edge]) {
                continue;
            }

            tm_triangle_edge_vertices(&mesh->triangles[tri_index], edge, &a, &b);
            source_segment = tm_find_snapshot_segment(snapshot, snapshot_count, a, b);
            if (source_segment != NULL) {
                original_index = source_segment->original_index;
            }

            status = tm_append_live_segment_unique(mesh, a, b, original_index, source_segment);
            if (status != TM_OK) {
                free(snapshot);
                return status;
            }
        }
    }

    free(snapshot);
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

static void tm_set_invalid_mesh_validation_detail(
    const TMMesh *mesh,
    const char *prefix
)
{
    const char *detail = tm_last_pslg_error_detail();
    TMValidationReport report;
    TMStatus validation_status;

    if (mesh == NULL || prefix == NULL) {
        return;
    }

    validation_status = tm_validate_mesh(mesh, 0, &report);
    tm_set_pslg_error_detail(
        "%s%s%s (validation status=%d adjacency=%zu orientation=%zu duplicates=%zu incident=%zu constrained=%zu)",
        prefix,
        (detail != NULL && detail[0] != '\0') ? ": " : "",
        (detail != NULL) ? detail : "",
        (int) validation_status,
        report.adjacency_errors,
        report.orientation_errors,
        report.duplicate_triangle_errors,
        report.incident_triangle_errors,
        report.constrained_edge_errors
    );
}

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

static TMStatus tm_grow_encroachment_candidates(
    TMEncroachmentCandidate **candidates,
    size_t *capacity,
    size_t min_capacity
)
{
    TMEncroachmentCandidate *grown;
    size_t new_capacity;

    if (*capacity >= min_capacity) {
        return TM_OK;
    }

    new_capacity = (*capacity == 0) ? 16 : *capacity;
    while (new_capacity < min_capacity) {
        if (new_capacity > ((size_t) -1) / 2) {
            return TM_ERR_ALLOC;
        }
        new_capacity *= 2;
    }

    grown = (TMEncroachmentCandidate *) realloc(*candidates, new_capacity * sizeof(*grown));
    if (grown == NULL) {
        return TM_ERR_ALLOC;
    }

    *candidates = grown;
    *capacity = new_capacity;
    return TM_OK;
}

static TMStatus tm_append_encroachment_candidate(
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity,
    size_t segment_index,
    int point_index
)
{
    TMStatus status = tm_grow_encroachment_candidates(candidates, capacity, *count + 1);

    if (status != TM_OK) {
        return status;
    }

    (*candidates)[*count].segment_index = segment_index;
    (*candidates)[*count].point_index = point_index;
    *count += 1;
    return TM_OK;
}

static TMStatus tm_grow_bad_triangle_flags(
    unsigned char **queued_flags,
    size_t *capacity,
    size_t min_capacity
)
{
    unsigned char *new_flags;
    size_t old_capacity;
    size_t new_capacity;

    if (*capacity >= min_capacity) {
        return TM_OK;
    }

    old_capacity = *capacity;
    new_capacity = (*capacity == 0) ? 32u : *capacity;
    while (new_capacity < min_capacity) {
        if (new_capacity > SIZE_MAX / 2u) {
            return TM_ERR_ALLOC;
        }
        new_capacity *= 2u;
    }

    new_flags = (unsigned char *) realloc(*queued_flags, new_capacity * sizeof(*new_flags));
    if (new_flags == NULL) {
        return TM_ERR_ALLOC;
    }

    memset(new_flags + old_capacity, 0, (new_capacity - old_capacity) * sizeof(*new_flags));
    *queued_flags = new_flags;
    *capacity = new_capacity;
    return TM_OK;
}

static TMStatus tm_grow_bad_triangle_candidates(
    TMBadTriangleCandidate **candidates,
    size_t *capacity,
    size_t min_capacity
)
{
    TMBadTriangleCandidate *new_candidates;
    size_t new_capacity;

    if (*capacity >= min_capacity) {
        return TM_OK;
    }

    new_capacity = (*capacity == 0) ? 32u : *capacity;
    while (new_capacity < min_capacity) {
        if (new_capacity > SIZE_MAX / 2u) {
            return TM_ERR_ALLOC;
        }
        new_capacity *= 2u;
    }

    new_candidates = (TMBadTriangleCandidate *) realloc(*candidates, new_capacity * sizeof(*new_candidates));
    if (new_candidates == NULL) {
        return TM_ERR_ALLOC;
    }

    *candidates = new_candidates;
    *capacity = new_capacity;
    return TM_OK;
}

static TMStatus tm_append_bad_triangle_candidate(
    TMBadTriangleCandidate **candidates,
    size_t *count,
    size_t *capacity,
    unsigned char **queued_flags,
    size_t *queued_flag_capacity,
    size_t triangle_index
)
{
    TMStatus status;

    status = tm_grow_bad_triangle_flags(queued_flags, queued_flag_capacity, triangle_index + 1);
    if (status != TM_OK) {
        return status;
    }
    if ((*queued_flags)[triangle_index]) {
        return TM_OK;
    }

    status = tm_grow_bad_triangle_candidates(candidates, capacity, *count + 1);

    if (status != TM_OK) {
        return status;
    }

    (*candidates)[*count].triangle_index = triangle_index;
    (*queued_flags)[triangle_index] = 1;
    *count += 1;
    return TM_OK;
}

static int tm_find_triangle_vertex_index(const TMTriangle *triangle, int vertex_index)
{
    int slot;

    for (slot = 0; slot < 3; ++slot) {
        if (triangle->v[slot] == vertex_index) {
            return slot;
        }
    }

    return -1;
}

static TMStatus tm_maybe_enqueue_bad_triangle_candidate(
    const TMMesh *mesh,
    size_t triangle_index,
    double beta,
    double max_area,
    TMBadTriangleCandidate **candidates,
    size_t *count,
    size_t *capacity,
    unsigned char **queued_flags,
    size_t *queued_flag_capacity
)
{
    if (!tm_evaluate_triangle_refinement_need(mesh, triangle_index, beta, max_area, NULL, NULL, NULL)) {
        return TM_OK;
    }

    return tm_append_bad_triangle_candidate(
        candidates,
        count,
        capacity,
        queued_flags,
        queued_flag_capacity,
        triangle_index
    );
}

static TMStatus tm_seed_bad_triangle_candidates(
    const TMMesh *mesh,
    double beta,
    double max_area,
    TMBadTriangleCandidate **candidates,
    size_t *count,
    size_t *capacity,
    unsigned char **queued_flags,
    size_t *queued_flag_capacity
)
{
    size_t triangle_index;

    for (triangle_index = 0; triangle_index < mesh->triangle_count; ++triangle_index) {
        TMStatus status = tm_maybe_enqueue_bad_triangle_candidate(
            mesh,
            triangle_index,
            beta,
            max_area,
            candidates,
            count,
            capacity,
            queued_flags,
            queued_flag_capacity
        );

        if (status != TM_OK) {
            return status;
        }
    }

    return TM_OK;
}

static TMStatus tm_enqueue_bad_triangles_along_point_fan(
    const TMMesh *mesh,
    int point_index,
    int triangle_index,
    int previous_triangle,
    int stop_triangle,
    double beta,
    double max_area,
    TMBadTriangleCandidate **candidates,
    size_t *count,
    size_t *capacity,
    unsigned char **queued_flags,
    size_t *queued_flag_capacity
)
{
    while (triangle_index >= 0) {
        const TMTriangle *triangle;
        int vertex_slot;
        int next_triangle;
        TMStatus status;

        if ((size_t) triangle_index >= mesh->triangle_count) {
            tm_set_pslg_error_detail(
                "enqueue bad triangles: point %d reached out-of-range triangle %d while traversing fan",
                point_index,
                triangle_index
            );
            return TM_ERR_INVALID_MESH;
        }

        triangle = &mesh->triangles[triangle_index];
        vertex_slot = tm_find_triangle_vertex_index(triangle, point_index);
        if (vertex_slot < 0) {
            tm_set_pslg_error_detail(
                "enqueue bad triangles: point %d is not a vertex of triangle %d during fan traversal",
                point_index,
                triangle_index
            );
            return TM_ERR_INVALID_MESH;
        }

        status = tm_maybe_enqueue_bad_triangle_candidate(
            mesh,
            (size_t) triangle_index,
            beta,
            max_area,
            candidates,
            count,
            capacity,
            queued_flags,
            queued_flag_capacity
        );
        if (status != TM_OK) {
            return status;
        }

        if (triangle->nbr[(vertex_slot + 1) % 3] == previous_triangle) {
            next_triangle = triangle->nbr[(vertex_slot + 2) % 3];
        } else if (triangle->nbr[(vertex_slot + 2) % 3] == previous_triangle) {
            next_triangle = triangle->nbr[(vertex_slot + 1) % 3];
        } else {
            tm_set_pslg_error_detail(
                "enqueue bad triangles: point %d triangle %d fan traversal lost predecessor %d",
                point_index,
                triangle_index,
                previous_triangle
            );
            return TM_ERR_INVALID_MESH;
        }

        if (next_triangle < 0 || next_triangle == stop_triangle) {
            return TM_OK;
        }

        previous_triangle = triangle_index;
        triangle_index = next_triangle;
    }

    return TM_OK;
}

static TMStatus tm_enqueue_bad_triangles_for_point(
    const TMMesh *mesh,
    int point_index,
    double beta,
    double max_area,
    TMBadTriangleCandidate **candidates,
    size_t *count,
    size_t *capacity,
    unsigned char **queued_flags,
    size_t *queued_flag_capacity
)
{
    int start_triangle;
    const TMTriangle *triangle;
    int vertex_slot;
    int first_neighbor;
    int second_neighbor;
    TMStatus status;

    if (point_index < 0 || (size_t) point_index >= mesh->point_count) {
        return TM_OK;
    }
    if (mesh->points[point_index].incident_triangle < 0 ||
        (size_t) mesh->points[point_index].incident_triangle >= mesh->triangle_count) {
        tm_set_pslg_error_detail(
            "enqueue bad triangles: point %d has invalid incident triangle %d",
            point_index,
            mesh->points[point_index].incident_triangle
        );
        return TM_ERR_INVALID_MESH;
    }

    start_triangle = mesh->points[point_index].incident_triangle;
    triangle = &mesh->triangles[start_triangle];
    vertex_slot = tm_find_triangle_vertex_index(triangle, point_index);
    if (vertex_slot < 0) {
        tm_set_pslg_error_detail(
            "enqueue bad triangles: point %d incident triangle %d does not contain the point",
            point_index,
            start_triangle
        );
        return TM_ERR_INVALID_MESH;
    }

    status = tm_maybe_enqueue_bad_triangle_candidate(
        mesh,
        (size_t) start_triangle,
        beta,
        max_area,
        candidates,
        count,
        capacity,
        queued_flags,
        queued_flag_capacity
    );
    if (status != TM_OK) {
        return status;
    }

    first_neighbor = triangle->nbr[(vertex_slot + 1) % 3];
    second_neighbor = triangle->nbr[(vertex_slot + 2) % 3];

    status = tm_enqueue_bad_triangles_along_point_fan(
        mesh,
        point_index,
        first_neighbor,
        start_triangle,
        start_triangle,
        beta,
        max_area,
        candidates,
        count,
        capacity,
        queued_flags,
        queued_flag_capacity
    );
    if (status != TM_OK) {
        return status;
    }

    if (second_neighbor == first_neighbor) {
        return TM_OK;
    }

    return tm_enqueue_bad_triangles_along_point_fan(
        mesh,
        point_index,
        second_neighbor,
        start_triangle,
        start_triangle,
        beta,
        max_area,
        candidates,
        count,
        capacity,
        queued_flags,
        queued_flag_capacity
    );
}

static int tm_point_encroaches_live_segment(const TMMesh *mesh, size_t segment_index, int point_index)
{
    const TMSegment *segment;
    int triangle_index;
    int edge;

    if (segment_index >= mesh->segment_count ||
        point_index < 0 ||
        (size_t) point_index >= mesh->point_count) {
        return 0;
    }

    segment = &mesh->segments[segment_index];
    if (!segment->live || segment->is_protected) {
        return 0;
    }
    if (point_index == segment->v[0] || point_index == segment->v[1]) {
        return 0;
    }
    if (mesh->points[point_index].incident_triangle < 0) {
        return 0;
    }
    if (!tm_mesh_has_edge(mesh, segment->v[0], segment->v[1], &triangle_index, &edge)) {
        return 0;
    }
    if (mesh->triangles[triangle_index].v[edge] != point_index) {
        int neighbor_index = mesh->triangles[triangle_index].nbr[edge];

        if (neighbor_index < 0) {
            return 0;
        }

        edge = tm_find_edge_in_triangle(&mesh->triangles[neighbor_index], segment->v[0], segment->v[1]);
        if (edge < 0 || mesh->triangles[neighbor_index].v[edge] != point_index) {
            return 0;
        }
    }

    return tm_point_encroaches_segment(
        mesh->points[point_index].xy,
        mesh->points[segment->v[0]].xy,
        mesh->points[segment->v[1]].xy,
        tm_segment_length_squared(mesh, segment->v[0], segment->v[1])
    );
}

static TMStatus tm_maybe_enqueue_encroached_segment_from_triangle(
    const TMMesh *mesh,
    int point_index,
    int triangle_index,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    const TMTriangle *triangle;
    int vertex_slot;
    int a;
    int b;
    size_t segment_index;

    if (triangle_index < 0 || (size_t) triangle_index >= mesh->triangle_count) {
        return TM_OK;
    }

    triangle = &mesh->triangles[triangle_index];
    vertex_slot = tm_find_triangle_vertex_index(triangle, point_index);
    if (vertex_slot < 0 || !triangle->constrained[vertex_slot]) {
        return TM_OK;
    }

    tm_triangle_edge_vertices(triangle, vertex_slot, &a, &b);
    if (!tm_find_live_segment_index(mesh, a, b, &segment_index)) {
        tm_set_pslg_error_detail(
            "enqueue encroached segment: constrained edge %d-%d opposite point %d is missing from live segment registry",
            a,
            b,
            point_index
        );
        return TM_ERR_INVALID_MESH;
    }
    if (!tm_point_encroaches_live_segment(mesh, segment_index, point_index)) {
        return TM_OK;
    }

    return tm_append_encroachment_candidate(candidates, count, capacity, segment_index, point_index);
}

static TMStatus tm_enqueue_segment_encroachments_along_point_fan(
    const TMMesh *mesh,
    int point_index,
    int triangle_index,
    int previous_triangle,
    int stop_triangle,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    while (triangle_index >= 0) {
        const TMTriangle *triangle;
        int vertex_slot;
        int next_triangle;
        TMStatus status;

        if ((size_t) triangle_index >= mesh->triangle_count) {
            tm_set_pslg_error_detail(
                "enqueue encroached segments: point %d reached out-of-range triangle %d while traversing fan",
                point_index,
                triangle_index
            );
            return TM_ERR_INVALID_MESH;
        }

        triangle = &mesh->triangles[triangle_index];
        vertex_slot = tm_find_triangle_vertex_index(triangle, point_index);
        if (vertex_slot < 0) {
            tm_set_pslg_error_detail(
                "enqueue encroached segments: point %d is not a vertex of triangle %d during fan traversal",
                point_index,
                triangle_index
            );
            return TM_ERR_INVALID_MESH;
        }

        status = tm_maybe_enqueue_encroached_segment_from_triangle(
            mesh,
            point_index,
            triangle_index,
            candidates,
            count,
            capacity
        );
        if (status != TM_OK) {
            return status;
        }

        if (triangle->nbr[(vertex_slot + 1) % 3] == previous_triangle) {
            next_triangle = triangle->nbr[(vertex_slot + 2) % 3];
        } else if (triangle->nbr[(vertex_slot + 2) % 3] == previous_triangle) {
            next_triangle = triangle->nbr[(vertex_slot + 1) % 3];
        } else {
            tm_set_pslg_error_detail(
                "enqueue encroached segments: point %d triangle %d fan traversal lost predecessor %d",
                point_index,
                triangle_index,
                previous_triangle
            );
            return TM_ERR_INVALID_MESH;
        }

        if (next_triangle < 0 || next_triangle == stop_triangle) {
            return TM_OK;
        }

        previous_triangle = triangle_index;
        triangle_index = next_triangle;
    }

    return TM_OK;
}

static TMStatus tm_enqueue_segment_encroachments_for_point(
    const TMMesh *mesh,
    int point_index,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    int start_triangle;
    const TMTriangle *triangle;
    int vertex_slot;
    int first_neighbor;
    int second_neighbor;
    TMStatus status;

    if (point_index < 0 || (size_t) point_index >= mesh->point_count) {
        return TM_OK;
    }
    if (mesh->points[point_index].incident_triangle < 0) {
        return TM_OK;
    }

    start_triangle = mesh->points[point_index].incident_triangle;
    triangle = &mesh->triangles[start_triangle];
    vertex_slot = tm_find_triangle_vertex_index(triangle, point_index);
    if (vertex_slot < 0) {
        tm_set_pslg_error_detail(
            "enqueue encroached segments: point %d incident triangle %d does not contain the point",
            point_index,
            start_triangle
        );
        return TM_ERR_INVALID_MESH;
    }

    status = tm_maybe_enqueue_encroached_segment_from_triangle(
        mesh,
        point_index,
        start_triangle,
        candidates,
        count,
        capacity
    );
    if (status != TM_OK) {
        return status;
    }

    first_neighbor = triangle->nbr[(vertex_slot + 1) % 3];
    second_neighbor = triangle->nbr[(vertex_slot + 2) % 3];

    status = tm_enqueue_segment_encroachments_along_point_fan(
        mesh,
        point_index,
        first_neighbor,
        start_triangle,
        start_triangle,
        candidates,
        count,
        capacity
    );
    if (status != TM_OK) {
        return status;
    }

    if (second_neighbor == first_neighbor) {
        return TM_OK;
    }

    return tm_enqueue_segment_encroachments_along_point_fan(
        mesh,
        point_index,
        second_neighbor,
        start_triangle,
        start_triangle,
        candidates,
        count,
        capacity
    );
}

static TMStatus tm_enqueue_point_encroachments_for_segment(
    const TMMesh *mesh,
    size_t segment_index,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    if (segment_index >= mesh->segment_count) {
        return TM_OK;
    }
    if (!mesh->segments[segment_index].live || mesh->segments[segment_index].is_protected) {
        return TM_OK;
    }

    {
        int triangle_index;
        int edge;

        if (!tm_mesh_has_edge(mesh, mesh->segments[segment_index].v[0], mesh->segments[segment_index].v[1], &triangle_index, &edge)) {
            return TM_OK;
        }

        if (tm_point_encroaches_live_segment(mesh, segment_index, mesh->triangles[triangle_index].v[edge])) {
            return tm_append_encroachment_candidate(
                candidates,
                count,
                capacity,
                segment_index,
                mesh->triangles[triangle_index].v[edge]
            );
        }

        if (mesh->triangles[triangle_index].nbr[edge] >= 0) {
            int neighbor_index = mesh->triangles[triangle_index].nbr[edge];
            int neighbor_edge = tm_find_edge_in_triangle(
                &mesh->triangles[neighbor_index],
                mesh->segments[segment_index].v[0],
                mesh->segments[segment_index].v[1]
            );

            if (neighbor_edge >= 0 &&
                tm_point_encroaches_live_segment(mesh, segment_index, mesh->triangles[neighbor_index].v[neighbor_edge])) {
                return tm_append_encroachment_candidate(
                    candidates,
                    count,
                    capacity,
                    segment_index,
                    mesh->triangles[neighbor_index].v[neighbor_edge]
                );
            }
        }
    }

    return TM_OK;
}

static TMStatus tm_seed_encroachment_candidates(
    const TMMesh *mesh,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    size_t segment_index;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        TMStatus status = tm_enqueue_point_encroachments_for_segment(
            mesh,
            segment_index,
            candidates,
            count,
            capacity
        );

        if (status != TM_OK) {
            return status;
        }
    }

    return TM_OK;
}

static int tm_pop_encroachment_candidate(
    const TMMesh *mesh,
    TMEncroachmentCandidate *candidates,
    size_t *count,
    size_t *out_segment_index,
    int *out_point_index
)
{
    while (*count > 0) {
        size_t candidate_index;
        size_t best_index = *count;
        double best_length_sq = 0.0;

        for (candidate_index = 0; candidate_index < *count; ++candidate_index) {
            TMEncroachmentCandidate candidate = candidates[candidate_index];
            double length_sq;

            if (!tm_point_encroaches_live_segment(mesh, candidate.segment_index, candidate.point_index)) {
                continue;
            }

            length_sq = tm_segment_length_squared(
                mesh,
                mesh->segments[candidate.segment_index].v[0],
                mesh->segments[candidate.segment_index].v[1]
            );
            if (best_index == *count || length_sq < best_length_sq) {
                best_index = candidate_index;
                best_length_sq = length_sq;
            }
        }

        if (best_index == *count) {
            *count = 0;
            return 0;
        }

        {
            TMEncroachmentCandidate candidate = candidates[best_index];

            candidates[best_index] = candidates[*count - 1];
            *count -= 1;
            if (!tm_point_encroaches_live_segment(mesh, candidate.segment_index, candidate.point_index)) {
                continue;
            }

            if (out_segment_index != NULL) {
                *out_segment_index = candidate.segment_index;
            }
            if (out_point_index != NULL) {
                *out_point_index = candidate.point_index;
            }
            return 1;
        }
    }

    return 0;
}

static TMStatus tm_enqueue_encroachments_after_segment_split(
    const TMMesh *mesh,
    int a,
    int b,
    int point_index,
    TMEncroachmentCandidate **candidates,
    size_t *count,
    size_t *capacity
)
{
    size_t left_segment_index;
    size_t right_segment_index;
    TMStatus status;

    if (point_index < 0) {
        return TM_OK;
    }

    status = tm_enqueue_segment_encroachments_for_point(
        mesh,
        point_index,
        candidates,
        count,
        capacity
    );
    if (status != TM_OK) {
        return status;
    }

    if (tm_find_live_segment_index(mesh, a, point_index, &left_segment_index)) {
        status = tm_enqueue_point_encroachments_for_segment(
            mesh,
            left_segment_index,
            candidates,
            count,
            capacity
        );
        if (status != TM_OK) {
            return status;
        }
    }

    if (tm_find_live_segment_index(mesh, point_index, b, &right_segment_index)) {
        status = tm_enqueue_point_encroachments_for_segment(
            mesh,
            right_segment_index,
            candidates,
            count,
            capacity
        );
        if (status != TM_OK) {
            return status;
        }
    }

    return TM_OK;
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

static size_t tm_collect_live_incident_segments(
    const TMMesh *mesh,
    int apex,
    size_t *out_segments,
    size_t out_capacity
)
{
    size_t segment_index;
    size_t count = 0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        if (!mesh->segments[segment_index].live) {
            continue;
        }
        if (mesh->segments[segment_index].v[0] == apex || mesh->segments[segment_index].v[1] == apex) {
            if (out_segments != NULL && count < out_capacity) {
                out_segments[count] = segment_index;
            }
            count += 1;
        }
    }

    return count;
}

static int tm_find_live_incident_segments(const TMMesh *mesh, int apex, size_t out_segments[2])
{
    return tm_collect_live_incident_segments(mesh, apex, out_segments, 2) == 2;
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

static int tm_find_smallest_acute_incident_pair(
    const TMMesh *mesh,
    int apex,
    const size_t *incident_segments,
    size_t incident_count,
    double protect_angle_deg,
    size_t out_segments[2],
    int *out_left,
    int *out_right,
    double *out_angle_deg
)
{
    size_t i;
    int found = 0;
    double best_angle = protect_angle_deg;

    for (i = 0; i < incident_count; ++i) {
        const TMSegment *first = &mesh->segments[incident_segments[i]];
        int first_other;
        size_t j;

        if (!first->live || first->is_protected) {
            continue;
        }

        first_other = tm_segment_other_endpoint(first, apex);
        for (j = i + 1; j < incident_count; ++j) {
            const TMSegment *second = &mesh->segments[incident_segments[j]];
            int second_other;
            double angle_deg;

            if (!second->live || second->is_protected) {
                continue;
            }

            second_other = tm_segment_other_endpoint(second, apex);
            if (first_other == second_other) {
                continue;
            }

            angle_deg = tm_corner_angle_deg(mesh, apex, first_other, second_other);
            if (!(angle_deg < protect_angle_deg * (1.0 - 1e-12))) {
                continue;
            }

            if (!found || angle_deg < best_angle) {
                found = 1;
                best_angle = angle_deg;
                out_segments[0] = incident_segments[i];
                out_segments[1] = incident_segments[j];
                if (out_left != NULL) {
                    *out_left = first_other;
                }
                if (out_right != NULL) {
                    *out_right = second_other;
                }
                if (out_angle_deg != NULL) {
                    *out_angle_deg = angle_deg;
                }
            }
        }
    }

    return found;
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

static void tm_copy_segment_flags(TMSegment *destination, const TMSegment *source)
{
    destination->is_protected = source->is_protected;
    destination->protected_apex = source->protected_apex;
}

static void tm_merge_segment_flags(TMSegment *destination, const TMSegment *source)
{
    if (source == NULL || !source->is_protected) {
        return;
    }

    destination->is_protected = 1;
    if (destination->protected_apex < 0) {
        destination->protected_apex = source->protected_apex;
    }
}

static TMStatus tm_append_live_segment_unique(
    TMMesh *mesh,
    int a,
    int b,
    int original_index,
    const TMSegment *source_flags
)
{
    size_t segment_index;
    TMStatus status;

    if (a == b) {
        return TM_OK;
    }

    if (tm_find_live_segment_index(mesh, a, b, &segment_index)) {
        tm_merge_segment_flags(&mesh->segments[segment_index], source_flags);
        return TM_OK;
    }

    segment_index = mesh->segment_count;
    status = tm_append_live_segment(mesh, a, b, original_index);
    if (status == TM_OK && source_flags != NULL) {
        tm_copy_segment_flags(&mesh->segments[segment_index], source_flags);
    }
    return status;
}

static TMStatus tm_split_segment_registry(TMMesh *mesh, size_t segment_index, int point_index)
{
    TMSegment segment = mesh->segments[segment_index];
    TMStatus status;

    mesh->segments[segment_index].live = 0;

    status = tm_append_live_segment_unique(
        mesh,
        segment.v[0],
        point_index,
        segment.original_index,
        &segment
    );
    if (status == TM_OK) {
        status = tm_append_live_segment_unique(
            mesh,
            point_index,
            segment.v[1],
            segment.original_index,
            &segment
        );
    }
    return status;
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
            TMStatus status = tm_append_live_segment_unique(
                target,
                a,
                point_index,
                segment->original_index,
                segment
            );
            if (status != TM_OK) {
                return status;
            }
            status = tm_append_live_segment_unique(
                target,
                point_index,
                b,
                segment->original_index,
                segment
            );
            if (status != TM_OK) {
                return status;
            }
            continue;
        }

        {
            TMStatus status = tm_append_live_segment_unique(
                target,
                segment->v[0],
                segment->v[1],
                segment->original_index,
                segment
            );
            if (status != TM_OK) {
                return status;
            }
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
    double bax = b[0] - a[0];
    double bay = b[1] - a[1];
    double cax = c[0] - a[0];
    double cay = c[1] - a[1];
    double balen_sq = bax * bax + bay * bay;
    double calen_sq = cax * cax + cay * cay;
    double det = 2.0 * (bax * cay - bay * cax);

    if (det == 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    out_point[0] = a[0] + (cay * balen_sq - bay * calen_sq) / det;
    out_point[1] = a[1] + (bax * calen_sq - cax * balen_sq) / det;
    return TM_OK;
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

static int tm_find_containing_triangle_bruteforce(
    const TMMesh *mesh,
    const double point[2],
    TMLocation *out_location
)
{
    size_t triangle_index;

    if (mesh == NULL || point == NULL || out_location == NULL) {
        return 0;
    }

    for (triangle_index = 0; triangle_index < mesh->triangle_count; ++triangle_index) {
        int on_edge = -1;

        if (!tm_point_inside_or_on_triangle(mesh, (int) triangle_index, point, &on_edge)) {
            continue;
        }

        out_location->triangle = (int) triangle_index;
        out_location->edge = on_edge;
        out_location->on_edge = (on_edge >= 0);
        return 1;
    }

    return 0;
}

static int tm_find_blocking_segment_for_point(
    const TMMesh *mesh,
    int start_triangle,
    const double point[2],
    size_t *out_segment_index
)
{
    TMLocation location;
    int on_edge = -1;
    TMStatus status;
    int a;
    int b;

    if (mesh == NULL || point == NULL || out_segment_index == NULL) {
        return 0;
    }

    status = tm_locate_point(mesh, point, start_triangle, &location);
    if (status != TM_OK) {
        location.triangle = -1;
        location.edge = -1;
    }
    if (status == TM_OK &&
        location.triangle >= 0 &&
        (size_t) location.triangle < mesh->triangle_count &&
        tm_point_inside_or_on_triangle(mesh, location.triangle, point, &on_edge)) {
        return 0;
    }
    if (location.triangle < 0 ||
        (size_t) location.triangle >= mesh->triangle_count ||
        location.edge < 0 ||
        location.edge > 2) {
        if (start_triangle >= 0 && (size_t) start_triangle < mesh->triangle_count) {
            const TMTriangle *triangle = &mesh->triangles[start_triangle];
            double centroid[2];
            size_t segment_index;
            int found = 0;
            double best_t = 0.0;
            double rx;
            double ry;

            centroid[0] = (
                mesh->points[triangle->v[0]].xy[0] +
                mesh->points[triangle->v[1]].xy[0] +
                mesh->points[triangle->v[2]].xy[0]
            ) / 3.0;
            centroid[1] = (
                mesh->points[triangle->v[0]].xy[1] +
                mesh->points[triangle->v[1]].xy[1] +
                mesh->points[triangle->v[2]].xy[1]
            ) / 3.0;
            rx = point[0] - centroid[0];
            ry = point[1] - centroid[1];

            for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
                const double *a;
                const double *b;
                double sx;
                double sy;
                double denom;
                double qpx;
                double qpy;
                double t;
                double u;

                if (!mesh->segments[segment_index].live) {
                    continue;
                }
                a = mesh->points[mesh->segments[segment_index].v[0]].xy;
                b = mesh->points[mesh->segments[segment_index].v[1]].xy;
                if (!tm_segments_cross_properly(
                        centroid,
                        point,
                        a,
                        b
                    )) {
                    continue;
                }
                sx = b[0] - a[0];
                sy = b[1] - a[1];
                denom = rx * sy - ry * sx;
                if (denom == 0.0) {
                    continue;
                }
                qpx = a[0] - centroid[0];
                qpy = a[1] - centroid[1];
                t = (qpx * sy - qpy * sx) / denom;
                u = (qpx * ry - qpy * rx) / denom;
                if (!(t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0)) {
                    continue;
                }
                if (!found || t < best_t) {
                    found = 1;
                    best_t = t;
                    *out_segment_index = segment_index;
                }
            }
            if (found) {
                return 1;
            }
        }
        return 0;
    }

    tm_triangle_edge_vertices(&mesh->triangles[location.triangle], location.edge, &a, &b);
    return tm_find_live_segment_index(mesh, a, b, out_segment_index);
}

static int tm_find_triangle_blocking_segment(
    const TMMesh *mesh,
    size_t triangle_index,
    const double point[2],
    size_t *out_segment_index
)
{
    const TMTriangle *triangle;
    int edge;
    int found = 0;
    double best_length_sq = 0.0;

    if (mesh == NULL || point == NULL || out_segment_index == NULL || triangle_index >= mesh->triangle_count) {
        return 0;
    }

    triangle = &mesh->triangles[triangle_index];
    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double side;
        int near_edge;
        size_t segment_index;
        double length_sq;

        if (!triangle->constrained[edge]) {
            continue;
        }

        tm_triangle_edge_vertices(triangle, edge, &a, &b);
        side = tm_orient_xy(mesh->points[a].xy, mesh->points[b].xy, point);
        near_edge = tm_point_nearly_on_segment_xy(mesh->points[a].xy, mesh->points[b].xy, point);
        if (side >= 0.0 || near_edge) {
            continue;
        }
        if (!tm_find_live_segment_index(mesh, a, b, &segment_index)) {
            continue;
        }

        length_sq = tm_segment_length_squared(mesh, a, b);
        if (!found || length_sq < best_length_sq) {
            found = 1;
            best_length_sq = length_sq;
            *out_segment_index = segment_index;
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

static int tm_find_protected_encroached_segment_for_point(
    const TMMesh *mesh,
    const double point[2],
    size_t *out_segment_index
)
{
    size_t segment_index;
    int found = 0;
    double best_length_sq = 0.0;

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        const TMSegment *segment = &mesh->segments[segment_index];
        double length_sq;

        if (!segment->live || !segment->is_protected) {
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

static int tm_evaluate_triangle_refinement_need(
    const TMMesh *mesh,
    size_t triangle_index,
    double beta,
    double max_area,
    double *out_ratio,
    double *out_area,
    int *out_quality_violation
)
{
    double ratio;
    int quality_exempt;

    if (triangle_index >= mesh->triangle_count) {
        return 0;
    }

    ratio = tm_triangle_radius_edge_ratio(mesh, (int) triangle_index);
    quality_exempt = tm_local_triangle_is_exempt_from_quality(mesh, (int) triangle_index);

    if (!quality_exempt && ratio > beta * (1.0 + 1e-12)) {
        if (out_ratio != NULL) {
            *out_ratio = ratio;
        }
        if (out_area != NULL) {
            *out_area = tm_triangle_area(mesh, (int) triangle_index);
        }
        if (out_quality_violation != NULL) {
            *out_quality_violation = 1;
        }
        return 1;
    }

    if (max_area > 0.0 && !quality_exempt) {
        double area = tm_triangle_area(mesh, (int) triangle_index);

        if (area > max_area * (1.0 + 1e-12)) {
            if (out_ratio != NULL) {
                *out_ratio = ratio;
            }
            if (out_area != NULL) {
                *out_area = area;
            }
            if (out_quality_violation != NULL) {
                *out_quality_violation = 0;
            }
            return 1;
        }
    }

    if (out_ratio != NULL) {
        *out_ratio = ratio;
    }
    if (out_area != NULL) {
        *out_area = tm_triangle_area(mesh, (int) triangle_index);
    }
    return 0;
}

static int tm_pop_bad_triangle_candidate(
    const TMMesh *mesh,
    TMBadTriangleCandidate *candidates,
    size_t *count,
    unsigned char *queued_flags,
    size_t queued_flag_capacity,
    double beta,
    double max_area,
    size_t *out_triangle_index,
    double *out_ratio,
    double *out_area,
    int *out_quality_violation
)
{
    while (*count != 0) {
        TMBadTriangleCandidate candidate = candidates[*count - 1];

        *count -= 1;
        if (candidate.triangle_index < queued_flag_capacity) {
            queued_flags[candidate.triangle_index] = 0;
        }
        if (tm_evaluate_triangle_refinement_need(
                mesh,
                candidate.triangle_index,
                beta,
                max_area,
                out_ratio,
                out_area,
                out_quality_violation
            )) {
            *out_triangle_index = candidate.triangle_index;
            return 1;
        }
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
    int use_existing_point = 0;
    int point_index;
    size_t rebuild_point_count;
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
    use_existing_point = blocker_index >= 0;
    point_index = use_existing_point ? blocker_index : (int) mesh->point_count;
    rebuild_point_count = mesh->point_count + (use_existing_point ? 0u : 1u);
    input_points = (TMPoint *) calloc(rebuild_point_count, sizeof(*input_points));
    if (input_points == NULL) {
        return TM_ERR_ALLOC;
    }

    memcpy(input_points, mesh->points, mesh->point_count * sizeof(*input_points));
    if (!use_existing_point) {
        input_points[point_index].xy[0] = split_point[0];
        input_points[point_index].xy[1] = split_point[1];
        input_points[point_index].original_index = point_index;
        input_points[point_index].kind = TM_VERTEX_SEGMENT_SPLIT;
        input_points[point_index].incident_triangle = -1;
        input_points[point_index].protection_apex = -1;
        input_points[point_index].protection_side = TM_PROTECTION_SIDE_NONE;
        input_points[point_index].protection_level = 0;
    }

    status = tm_build_mesh(input_points, rebuild_point_count, &rebuilt);
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
        for (i = 0; i < rebuild_point_count; ++i) {
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
            tm_verbose_log(
                options,
                "recover segment %d-%d: edge flip failed with status %d; split and retry",
                a,
                b,
                (int) status
            );
            break;
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
                if (!tm_edge_is_flippable(mesh, (int) tri_index, edge)) {
                    continue;
                }

                if (!tm_edge_is_locally_delaunay(mesh, (int) tri_index, edge)) {
                    TMStatus status;
                    const char *detail;

                    tm_verbose_log(options, "restore CDT: flip edge in triangle %zu edge %d", tri_index, edge);
                    status = tm_flip_edge(mesh, (int) tri_index, edge);
                    if (status != TM_OK) {
                        detail = tm_last_pslg_error_detail();
                        {
                            TMValidationReport report;
                            TMStatus validation_status = tm_validate_mesh(mesh, 0, &report);

                            tm_set_pslg_error_detail(
                                "restore constrained Delaunay failed at triangle %zu edge %d%s%s (validation status=%d adjacency=%zu orientation=%zu duplicates=%zu incident=%zu constrained=%zu)",
                                tri_index,
                                edge,
                                (detail != NULL && detail[0] != '\0') ? ": " : "",
                                (detail != NULL) ? detail : "",
                                (int) validation_status,
                                report.adjacency_errors,
                                report.orientation_errors,
                                report.duplicate_triangle_errors,
                                report.incident_triangle_errors,
                                report.constrained_edge_errors
                            );
                        }
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
    const char *detail;
    double length_sq;
    double tolerance;

    (void) options;

    if (!tm_mesh_has_edge(mesh, segment->v[0], segment->v[1], &triangle_index, &edge)) {
        tm_set_pslg_error_detail(
            "split segment: live mesh edge %d-%d is missing for segment index %zu",
            segment->v[0],
            segment->v[1],
            segment_index
        );
        return TM_ERR_INVALID_MESH;
    }
    length_sq = tm_segment_length_squared(mesh, segment->v[0], segment->v[1]);
    tolerance = tm_segment_split_tolerance_xy(
        mesh->points[segment->v[0]].xy,
        mesh->points[segment->v[1]].xy,
        point
    );
    status = tm_insert_point_on_edge(mesh, triangle_index, edge, point, kind, &point_index);
    if (status != TM_OK) {
        if (length_sq <= 4096.0 * tolerance * tolerance) {
            mesh->segments[segment_index].is_protected = 1;
            if (out_point_index != NULL) {
                *out_point_index = -1;
            }
            tm_verbose_log(
                options,
                "split segment: protecting numerically unresolved segment %d-%d after edge split failure",
                segment->v[0],
                segment->v[1]
            );
            return TM_OK;
        }
        detail = tm_last_pslg_error_detail();
        tm_set_pslg_error_detail(
            "split segment %d-%d failed while inserting point (%.12g, %.12g), len=%.12g%s%s",
            segment->v[0],
            segment->v[1],
            point[0],
            point[1],
            sqrt(tm_segment_length_squared(mesh, segment->v[0], segment->v[1])),
            (detail != NULL && detail[0] != '\0') ? ": " : "",
            (detail != NULL) ? detail : ""
        );
        return status;
    }

    status = tm_split_segment_registry(mesh, segment_index, point_index);
    if (status != TM_OK) {
        tm_set_pslg_error_detail(
            "split segment: registry update failed for %d-%d at point %d",
            segment->v[0],
            segment->v[1],
            point_index
        );
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
    double length_sq;
    double tolerance;

    midpoint[0] = 0.5 * (mesh->points[segment->v[0]].xy[0] + mesh->points[segment->v[1]].xy[0]);
    midpoint[1] = 0.5 * (mesh->points[segment->v[0]].xy[1] + mesh->points[segment->v[1]].xy[1]);
    length_sq = tm_segment_length_squared(mesh, segment->v[0], segment->v[1]);
    tolerance = tm_segment_split_tolerance_xy(
        mesh->points[segment->v[0]].xy,
        mesh->points[segment->v[1]].xy,
        midpoint
    );
    if (length_sq <= 4.0 * tolerance * tolerance) {
        mesh->segments[segment_index].is_protected = 1;
        if (out_point_index != NULL) {
            *out_point_index = -1;
        }
        tm_verbose_log(
            options,
            "split segment midpoint: protecting numerically unresolved segment %d-%d",
            segment->v[0],
            segment->v[1]
        );
        return TM_OK;
    }
    {
        TMStatus status = tm_split_segment_with_point(
            mesh,
            segment_index,
            midpoint,
            TM_VERTEX_SEGMENT_SPLIT,
            options,
            out_point_index
        );

        if (status == TM_OK) {
            return TM_OK;
        }
        if (length_sq <= 4096.0 * tolerance * tolerance) {
            mesh->segments[segment_index].is_protected = 1;
            if (out_point_index != NULL) {
                *out_point_index = -1;
            }
            tm_verbose_log(
                options,
                "split segment midpoint: protecting numerically unresolved segment %d-%d after midpoint split failure",
                segment->v[0],
                segment->v[1]
            );
            return TM_OK;
        }
        return status;
    }
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
    int inner = apex;
    size_t level;

    if (segment_length <= 0.0) {
        return TM_ERR_INVALID_MESH;
    }

    for (level = 1; level <= levels; ++level) {
        double distance = base_length * (double) (1ULL << (level - 1));
        double t = distance / segment_length;
        double point[2];
        int point_index;
        TMStatus status;

        if (!(distance < segment_length * (1.0 - 1e-12))) {
            return TM_ERR_INVALID_MESH;
        }

        point[0] = mesh->points[apex].xy[0] + (mesh->points[outer].xy[0] - mesh->points[apex].xy[0]) * t;
        point[1] = mesh->points[apex].xy[1] + (mesh->points[outer].xy[1] - mesh->points[apex].xy[1]) * t;

        status = tm_split_segment_with_point(mesh, segment_index, point, TM_VERTEX_SEGMENT_SPLIT, options, &point_index);
        if (status != TM_OK) {
            const char *detail = tm_last_pslg_error_detail();

            tm_set_pslg_error_detail(
                "acute shell protection failed at apex %d level %zu on side %d%s%s",
                apex,
                level,
                (int) side,
                (detail != NULL && detail[0] != '\0') ? ": " : "",
                (detail != NULL) ? detail : ""
            );
            return status;
        }

        tm_mark_shell_point(mesh, point_index, apex, side, (unsigned int) level);
        status = tm_mark_live_segment_protected_with_apex(mesh, inner, point_index, apex);
        if (status != TM_OK) {
            tm_set_pslg_error_detail(
                "acute shell protection failed to protect segment %d-%d for apex %d",
                inner,
                point_index,
                apex
            );
            return status;
        }

        inner = point_index;
        if (level < levels && !tm_find_live_segment_index(mesh, inner, outer, &segment_index)) {
            tm_set_pslg_error_detail(
                "acute shell protection lost live segment %d-%d for apex %d",
                inner,
                outer,
                apex
            );
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

    status = tm_build_shell_chain(
        mesh,
        apex,
        right,
        incident_segments[1],
        TM_PROTECTION_SIDE_RIGHT,
        base_length,
        levels,
        options
    );
    if (status != TM_OK) {
        return status;
    }
    return TM_OK;
}

static TMStatus tm_apply_acute_corner_protection(TMMesh *mesh, const TMBuildOptions *options)
{
    double protect_angle_deg = tm_protection_angle_deg(options);
    size_t protection_level_limit = tm_protection_level_limit(options);
    size_t *incident_segments = NULL;
    int apex;
    TMStatus final_status = TM_OK;

    if (!tm_protect_acute_corners(options)) {
        return TM_OK;
    }

    if (mesh->segment_count != 0) {
        incident_segments = (size_t *) malloc(mesh->segment_count * sizeof(*incident_segments));
        if (incident_segments == NULL) {
            return TM_ERR_ALLOC;
        }
    }

    for (apex = 0; (size_t) apex < mesh->point_count; ++apex) {
        size_t incident_count;

        if (mesh->points[apex].kind != TM_VERTEX_INPUT) {
            continue;
        }

        while (1) {
            size_t protected_pair[2];
            int left;
            int right;
            double angle_deg;
            size_t levels = 1;
            TMStatus status;

            incident_count = tm_collect_live_incident_segments(mesh, apex, incident_segments, mesh->segment_count);
            if (incident_count < 2) {
                break;
            }
            if (!tm_find_smallest_acute_incident_pair(
                    mesh,
                    apex,
                    incident_segments,
                    incident_count,
                    protect_angle_deg,
                    protected_pair,
                    &left,
                    &right,
                    &angle_deg
                )) {
                break;
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
                status = tm_apply_simple_acute_corner_protection(mesh, protected_pair, apex, levels, options);
            } else {
                status = tm_apply_shell_acute_corner_protection(mesh, protected_pair, apex, left, right, levels, options);
            }
            if (status != TM_OK) {
                const char *detail = tm_last_pslg_error_detail();

                tm_set_pslg_error_detail(
                    "acute protection failed at apex %d%s%s",
                    apex,
                    (detail != NULL && detail[0] != '\0') ? ": " : "",
                    (detail != NULL) ? detail : ""
                );
                final_status = status;
                goto cleanup;
            }
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

cleanup:
    free(incident_segments);
    return final_status;
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
        if (tm_find_containing_triangle_bruteforce(mesh, point, &location)) {
            status = TM_OK;
        }
    }
    if (status != TM_OK) {
        const char *detail = tm_last_pslg_error_detail();
        size_t blocking_segment_index = 0;
        int blocking_found = tm_find_triangle_blocking_segment(mesh, (size_t) start_triangle, point, &blocking_segment_index);

        tm_set_pslg_error_detail(
            "insert triangle split point: locate failed from triangle %d at (%.12g, %.12g); vertices=%d,%d,%d constrained=%d,%d,%d blocking=%d%s%s",
            start_triangle,
            point[0],
            point[1],
            mesh->triangles[start_triangle].v[0],
            mesh->triangles[start_triangle].v[1],
            mesh->triangles[start_triangle].v[2],
            mesh->triangles[start_triangle].constrained[0],
            mesh->triangles[start_triangle].constrained[1],
            mesh->triangles[start_triangle].constrained[2],
            blocking_found ? (int) blocking_segment_index : -1,
            (detail != NULL && detail[0] != '\0') ? ": " : "",
            (detail != NULL) ? detail : ""
        );
        return status;
    }

    if (!tm_point_inside_or_on_triangle(mesh, location.triangle, point, &on_edge)) {
        if (!tm_find_containing_triangle_bruteforce(mesh, point, &location) ||
            !tm_point_inside_or_on_triangle(mesh, location.triangle, point, &on_edge)) {
            tm_set_pslg_error_detail(
                "insert triangle split point: located triangle %d does not contain point (%.12g, %.12g) from start triangle %d",
                location.triangle,
                point[0],
                point[1],
                start_triangle
            );
            return TM_ERR_INVALID_MESH;
        }
    }

    if (on_edge >= 0) {
        int a;
        int b;
        size_t segment_index;

        if (mesh->triangles[location.triangle].constrained[on_edge]) {
            tm_triangle_edge_vertices(&mesh->triangles[location.triangle], on_edge, &a, &b);
            if (!tm_find_live_segment_index(mesh, a, b, &segment_index)) {
                tm_set_pslg_error_detail(
                    "insert triangle split point: constrained edge %d-%d is missing from live segment registry",
                    a,
                    b
                );
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

    return TM_OK;
}

static TMStatus tm_refine_quality_mesh(TMMesh *mesh, const TMBuildOptions *options)
{
    double min_angle_deg = tm_quality_min_angle_deg(options);
    double max_area = tm_area_limit(options);
    double beta = tm_quality_beta(min_angle_deg);
    TMEncroachmentCandidate *encroachment_candidates = NULL;
    TMBadTriangleCandidate *bad_triangle_candidates = NULL;
    unsigned char *bad_triangle_queued = NULL;
    size_t encroachment_count = 0;
    size_t encroachment_capacity = 0;
    size_t bad_triangle_count = 0;
    size_t bad_triangle_capacity = 0;
    size_t bad_triangle_queued_capacity = 0;
    size_t step_limit = tm_refinement_step_limit(options, mesh);
    size_t steps = 0;
    size_t segment_splits = 0;
    size_t triangle_splits = 0;
    double seed_time = 0.0;
    double pop_time = 0.0;
    double bad_triangle_time = 0.0;
    double split_time = 0.0;
    double enqueue_time = 0.0;
    double circumcenter_time = 0.0;
    double triangle_insert_time = 0.0;
    TMStatus status;

    tm_clear_pslg_error_detail();

    if (max_area > 0.0) {
        tm_verbose_log(
            options,
            "quality refinement: target min angle %.2f deg (beta=%.6f), max area=%.6f, step limit=%zu, mode=circumcenters",
            min_angle_deg,
            beta,
            max_area,
            step_limit
        );
    } else {
        tm_verbose_log(
            options,
            "quality refinement: target min angle %.2f deg (beta=%.6f), max area=off, step limit=%zu, mode=circumcenters",
            min_angle_deg,
            beta,
            step_limit
        );
    }

    {
        clock_t start = clock();

        status = tm_seed_encroachment_candidates(
            mesh,
            &encroachment_candidates,
            &encroachment_count,
            &encroachment_capacity
        );
        if (status == TM_OK) {
            status = tm_seed_bad_triangle_candidates(
                mesh,
                beta,
                max_area,
                &bad_triangle_candidates,
                &bad_triangle_count,
                &bad_triangle_capacity,
                &bad_triangle_queued,
                &bad_triangle_queued_capacity
            );
        }
        seed_time += tm_elapsed_seconds(start, clock());
    }
    if (status != TM_OK) {
        free(bad_triangle_queued);
        free(bad_triangle_candidates);
        free(encroachment_candidates);
        return status;
    }

    while (steps < step_limit) {
        size_t segment_index;
        int encroaching_point = -1;
        size_t triangle_index;
        double ratio;
        double area;
        int quality_violation = 0;

        {
            clock_t start = clock();
            int have_encroached_segment = tm_pop_encroachment_candidate(
                mesh,
                encroachment_candidates,
                &encroachment_count,
                &segment_index,
                &encroaching_point
            );

            pop_time += tm_elapsed_seconds(start, clock());

            if (have_encroached_segment) {
                const TMSegment *segment = &mesh->segments[segment_index];
                int a = segment->v[0];
                int b = segment->v[1];
                int split_point_index = -1;

                tm_verbose_log(
                    options,
                    "refine step %zu: split encroached segment %d-%d due to point %d",
                    steps + 1,
                    segment->v[0],
                    segment->v[1],
                    encroaching_point
                );
                start = clock();
                status = tm_split_segment_midpoint(mesh, segment_index, options, &split_point_index);
                split_time += tm_elapsed_seconds(start, clock());
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while splitting encroached segment"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                start = clock();
                status = tm_enqueue_encroachments_after_segment_split(
                    mesh,
                    a,
                    b,
                    split_point_index,
                    &encroachment_candidates,
                    &encroachment_count,
                    &encroachment_capacity
                );
                enqueue_time += tm_elapsed_seconds(start, clock());
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while enqueuing encroachments after segment split"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                start = clock();
                status = tm_enqueue_bad_triangles_for_point(
                    mesh,
                    split_point_index,
                    beta,
                    max_area,
                    &bad_triangle_candidates,
                    &bad_triangle_count,
                    &bad_triangle_capacity,
                    &bad_triangle_queued,
                    &bad_triangle_queued_capacity
                );
                enqueue_time += tm_elapsed_seconds(start, clock());
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while enqueuing bad triangles after segment split"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                segment_splits += 1;
                steps += 1;
                continue;
            }
        }

        {
            clock_t start = clock();
            int found_bad_triangle = tm_pop_bad_triangle_candidate(
                mesh,
                bad_triangle_candidates,
                &bad_triangle_count,
                bad_triangle_queued,
                bad_triangle_queued_capacity,
                beta,
                max_area,
                &triangle_index,
                &ratio,
                &area,
                &quality_violation
            );

            bad_triangle_time += tm_elapsed_seconds(start, clock());

            if (!found_bad_triangle) {
                tm_verbose_log(
                    options,
                    "quality refinement complete: %zu steps, %zu segment splits, %zu triangle splits",
                    steps,
                    segment_splits,
                    triangle_splits
                );
                if (options != NULL && options->verbose) {
                    tm_verbose_log(
                        options,
                        "quality refinement timings: seed=%.3fs pop=%.3fs bad=%.3fs split=%.3fs enqueue=%.3fs circumcenter=%.3fs triangle_insert=%.3fs",
                        seed_time,
                        pop_time,
                        bad_triangle_time,
                        split_time,
                        enqueue_time,
                        circumcenter_time,
                        triangle_insert_time
                    );
                }
                free(bad_triangle_queued);
                free(bad_triangle_candidates);
                free(encroachment_candidates);
                return TM_OK;
            }
        }

        {
            double circumcenter[2];
            size_t circumcenter_segment_index = 0;
            size_t split_segment_index = 0;
            int circumcenter_encroaches = 0;
            int circumcenter_encroaches_protected = 0;
            int split_segment = 0;
            int point_index = -1;

            {
                clock_t start = clock();

                status = tm_triangle_circumcenter(mesh, (int) triangle_index, circumcenter);
                circumcenter_time += tm_elapsed_seconds(start, clock());
            }
            if (status != TM_OK) {
                free(bad_triangle_queued);
                free(bad_triangle_candidates);
                free(encroachment_candidates);
                return status;
            }

            circumcenter_encroaches = tm_find_encroached_segment_for_point(mesh, circumcenter, &circumcenter_segment_index);
            if (circumcenter_encroaches) {
                split_segment = 1;
                split_segment_index = circumcenter_segment_index;
            }

            if (!split_segment) {
                circumcenter_encroaches_protected = tm_find_protected_encroached_segment_for_point(
                    mesh,
                    circumcenter,
                    &circumcenter_segment_index
                );
                if (circumcenter_encroaches_protected) {
                    tm_verbose_log(
                        options,
                        "refine step %zu: skip bad triangle %zu because its circumcenter encroaches protected segment %d-%d",
                        steps + 1,
                        triangle_index,
                        mesh->segments[circumcenter_segment_index].v[0],
                        mesh->segments[circumcenter_segment_index].v[1]
                    );
                    continue;
                }
            }

            if (!split_segment &&
                tm_find_triangle_blocking_segment(mesh, triangle_index, circumcenter, &circumcenter_segment_index)) {
                if (mesh->segments[circumcenter_segment_index].is_protected) {
                    tm_verbose_log(
                        options,
                        "refine step %zu: skip bad triangle %zu because its circumcenter crosses protected segment %d-%d",
                        steps + 1,
                        triangle_index,
                        mesh->segments[circumcenter_segment_index].v[0],
                        mesh->segments[circumcenter_segment_index].v[1]
                    );
                    continue;
                }
                split_segment = 1;
                split_segment_index = circumcenter_segment_index;
            }

            if (!split_segment &&
                tm_find_blocking_segment_for_point(mesh, (int) triangle_index, circumcenter, &circumcenter_segment_index)) {
                if (mesh->segments[circumcenter_segment_index].is_protected) {
                    tm_verbose_log(
                        options,
                        "refine step %zu: skip bad triangle %zu because its circumcenter is blocked by protected segment %d-%d",
                        steps + 1,
                        triangle_index,
                        mesh->segments[circumcenter_segment_index].v[0],
                        mesh->segments[circumcenter_segment_index].v[1]
                    );
                    continue;
                }
                split_segment = 1;
                split_segment_index = circumcenter_segment_index;
            }

            if (split_segment) {
                const TMSegment *segment = &mesh->segments[split_segment_index];
                int a = segment->v[0];
                int b = segment->v[1];
                int split_point_index = -1;

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
                {
                    clock_t start = clock();

                    status = tm_split_segment_midpoint(mesh, split_segment_index, options, &split_point_index);
                    split_time += tm_elapsed_seconds(start, clock());
                }
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while splitting segment for bad triangle"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                {
                    clock_t start = clock();

                    status = tm_enqueue_encroachments_after_segment_split(
                        mesh,
                        a,
                        b,
                        split_point_index,
                        &encroachment_candidates,
                        &encroachment_count,
                        &encroachment_capacity
                    );
                    enqueue_time += tm_elapsed_seconds(start, clock());
                }
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while enqueuing encroachments after bad-triangle split"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                {
                    clock_t start = clock();

                    status = tm_enqueue_bad_triangles_for_point(
                        mesh,
                        split_point_index,
                        beta,
                        max_area,
                        &bad_triangle_candidates,
                        &bad_triangle_count,
                        &bad_triangle_capacity,
                        &bad_triangle_queued,
                        &bad_triangle_queued_capacity
                    );
                    if (status == TM_OK) {
                        status = tm_maybe_enqueue_bad_triangle_candidate(
                            mesh,
                            triangle_index,
                            beta,
                            max_area,
                            &bad_triangle_candidates,
                            &bad_triangle_count,
                            &bad_triangle_capacity,
                            &bad_triangle_queued,
                            &bad_triangle_queued_capacity
                        );
                    }
                    enqueue_time += tm_elapsed_seconds(start, clock());
                }
                if (status != TM_OK) {
                    if (status == TM_ERR_INVALID_MESH) {
                        tm_set_invalid_mesh_validation_detail(
                            mesh,
                            "quality refinement failed while enqueuing bad triangles after bad-triangle split"
                        );
                    }
                    free(bad_triangle_queued);
                    free(bad_triangle_candidates);
                    free(encroachment_candidates);
                    return status;
                }
                segment_splits += 1;
                steps += 1;
                continue;
            }

            tm_verbose_log(
                options,
                quality_violation ?
                    "refine step %zu: split bad triangle %zu at circumcenter (R/s=%.6f, area=%.6f)" :
                    "refine step %zu: split oversized triangle %zu at circumcenter (area=%.6f, R/s=%.6f)",
                steps + 1,
                triangle_index,
                quality_violation ? ratio : area,
                quality_violation ? area : ratio
            );
            {
                clock_t start = clock();
                size_t point_count_before = mesh->point_count;

                status = tm_insert_triangle_split_point(
                    mesh,
                    (int) triangle_index,
                    circumcenter,
                    options,
                    &point_index
                );
                triangle_insert_time += tm_elapsed_seconds(start, clock());
                if (status == TM_OK && mesh->point_count == point_count_before) {
                    tm_verbose_log(
                        options,
                        "refine step %zu: skip bad triangle %zu because its circumcenter coincides with an existing vertex",
                        steps + 1,
                        triangle_index
                    );
                    continue;
                }
            }
            if (status != TM_OK) {
                if (status == TM_ERR_INVALID_MESH) {
                    tm_set_invalid_mesh_validation_detail(
                        mesh,
                        "quality refinement failed while inserting triangle split point"
                    );
                }
                free(bad_triangle_queued);
                free(bad_triangle_candidates);
                free(encroachment_candidates);
                return status;
            }
            {
                clock_t start = clock();

                status = tm_enqueue_segment_encroachments_for_point(
                    mesh,
                    point_index,
                    &encroachment_candidates,
                    &encroachment_count,
                    &encroachment_capacity
                );
                if (status == TM_OK) {
                    status = tm_enqueue_bad_triangles_for_point(
                        mesh,
                        point_index,
                        beta,
                        max_area,
                        &bad_triangle_candidates,
                        &bad_triangle_count,
                        &bad_triangle_capacity,
                        &bad_triangle_queued,
                        &bad_triangle_queued_capacity
                    );
                }
                enqueue_time += tm_elapsed_seconds(start, clock());
            }
            if (status != TM_OK) {
                if (status == TM_ERR_INVALID_MESH) {
                    tm_set_invalid_mesh_validation_detail(
                        mesh,
                        "quality refinement failed while enqueuing after triangle split"
                    );
                }
                free(bad_triangle_queued);
                free(bad_triangle_candidates);
                free(encroachment_candidates);
                return status;
            }
            if (point_index >= 0 && mesh->points[point_index].kind == TM_VERTEX_SEGMENT_SPLIT) {
                size_t incident_segments[2];

                if (tm_find_live_incident_segments(mesh, point_index, incident_segments)) {
                    {
                        clock_t start = clock();

                        status = tm_enqueue_point_encroachments_for_segment(
                            mesh,
                            incident_segments[0],
                            &encroachment_candidates,
                            &encroachment_count,
                            &encroachment_capacity
                        );
                        enqueue_time += tm_elapsed_seconds(start, clock());
                    }
                    if (status != TM_OK) {
                        if (status == TM_ERR_INVALID_MESH) {
                            tm_set_invalid_mesh_validation_detail(
                                mesh,
                                "quality refinement failed while enqueuing point encroachments for first incident segment"
                            );
                        }
                        free(bad_triangle_queued);
                        free(bad_triangle_candidates);
                        free(encroachment_candidates);
                        return status;
                    }
                    {
                        clock_t start = clock();

                        status = tm_enqueue_point_encroachments_for_segment(
                            mesh,
                            incident_segments[1],
                            &encroachment_candidates,
                            &encroachment_count,
                            &encroachment_capacity
                        );
                        enqueue_time += tm_elapsed_seconds(start, clock());
                    }
                    if (status != TM_OK) {
                        if (status == TM_ERR_INVALID_MESH) {
                            tm_set_invalid_mesh_validation_detail(
                                mesh,
                                "quality refinement failed while enqueuing point encroachments for second incident segment"
                            );
                        }
                        free(bad_triangle_queued);
                        free(bad_triangle_candidates);
                        free(encroachment_candidates);
                        return status;
                    }
                }
            }
            triangle_splits += 1;
            steps += 1;
        }
    }

    free(bad_triangle_queued);
    free(bad_triangle_candidates);
    free(encroachment_candidates);
    if (options != NULL && options->verbose) {
        tm_verbose_log(
            options,
            "quality refinement timings: seed=%.3fs pop=%.3fs bad=%.3fs split=%.3fs enqueue=%.3fs circumcenter=%.3fs triangle_insert=%.3fs",
            seed_time,
            pop_time,
            bad_triangle_time,
            split_time,
            enqueue_time,
            circumcenter_time,
            triangle_insert_time
        );
    }
    tm_set_pslg_error_detail(
        "quality refinement reached step limit after %zu steps (%zu segment splits, %zu triangle splits)",
        steps,
        segment_splits,
        triangle_splits
    );
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

    return tm_sync_live_segments_from_constraints(mesh);
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
        int on_edge = -1;

        if (status != TM_OK ||
            location.triangle < 0 ||
            (size_t) location.triangle >= mesh->triangle_count ||
            !tm_point_inside_or_on_triangle(mesh, location.triangle, regions[region_index].xy, &on_edge)) {
            if (tm_find_containing_triangle_bruteforce(mesh, regions[region_index].xy, &location)) {
                status = TM_OK;
            }
        }

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

    return tm_sync_live_segments_from_constraints(mesh);
}

TMStatus tm_build_pslg_mesh(const TMPSLG *pslg, const TMBuildOptions *options, TMMesh *out_mesh)
{
    TMPSLG working_pslg;
    size_t segment_index;
    TMStatus status;

    if (pslg == NULL || out_mesh == NULL) {
        return TM_ERR_INTERNAL;
    }

    tm_clear_pslg_error_detail();
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
    int *pre_refine_markers = NULL;
    size_t segment_index;
    TMStatus status;

    if (pslg == NULL || out_mesh == NULL || out_triangle_markers == NULL) {
        return TM_ERR_INTERNAL;
    }
    if (regions == NULL || region_count == 0) {
        tm_set_pslg_error_detail("coverage meshing requires at least one region seed");
        return TM_ERR_INVALID_PSLG;
    }

    tm_clear_pslg_error_detail();
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

    tm_verbose_log(options, "coverage segment recovery complete");

    if (tm_refinement_enabled(options)) {
        status = tm_assign_coverage_regions(out_mesh, regions, region_count, &pre_refine_markers);
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }

        status = tm_compact_marked_triangles(out_mesh, &pre_refine_markers);
        free(pre_refine_markers);
        pre_refine_markers = NULL;
        if (status != TM_OK) {
            tm_free_pslg(&working_pslg);
            tm_free_mesh(out_mesh);
            return status;
        }

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
