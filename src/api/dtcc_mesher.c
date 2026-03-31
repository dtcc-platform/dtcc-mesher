#include "dtcc_mesher/dtcc_mesher.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cdt.h"
#include "io_pslg.h"
#include "mesh.h"

#include "dtcc_mesher_api_internal.h"

static void tm_api_init_point(TMPoint *point, double x, double y, int original_index, TMVertexKind kind)
{
    point->xy[0] = x;
    point->xy[1] = y;
    point->original_index = original_index;
    point->kind = kind;
    point->incident_triangle = -1;
    point->protection_apex = -1;
    point->protection_side = 0;
    point->protection_level = 0;
}

static void tm_api_init_segment(TMSegment *segment, uint32_t a, uint32_t b, int original_index)
{
    segment->v[0] = (int) a;
    segment->v[1] = (int) b;
    segment->original_index = original_index;
    segment->live = 1;
    segment->is_protected = 0;
    segment->protected_apex = -1;
}

static tm_status tm_validate_domain_shape(const tm_domain *domain, tm_error *out_error)
{
    if (domain == NULL) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "domain must not be null");
        return TM_STATUS_INVALID_ARGUMENT;
    }
    if (domain->num_points != 0 && domain->points == NULL) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "domain points must not be null when num_points > 0");
        return TM_STATUS_INVALID_ARGUMENT;
    }
    if (domain->num_segments != 0 && domain->segments == NULL) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "domain segments must not be null when num_segments > 0");
        return TM_STATUS_INVALID_ARGUMENT;
    }
    if (domain->num_holes != 0 && domain->holes == NULL) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "domain holes must not be null when num_holes > 0");
        return TM_STATUS_INVALID_ARGUMENT;
    }
    if (domain->num_segments == 0 && domain->num_holes != 0) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "holes require a segmented PSLG domain");
        return TM_STATUS_INVALID_ARGUMENT;
    }

    return TM_STATUS_OK;
}

static tm_status tm_copy_points_to_internal(const tm_domain *domain, TMPoint **out_points, tm_error *out_error)
{
    TMPoint *points = NULL;
    size_t i;

    *out_points = NULL;
    if (domain->num_points == 0) {
        tm_api_set_error(out_error, TM_STATUS_GEOMETRY, "at least three points are required");
        return TM_STATUS_GEOMETRY;
    }

    points = (TMPoint *) calloc(domain->num_points, sizeof(*points));
    if (points == NULL) {
        tm_api_set_error(out_error, TM_STATUS_NO_MEMORY, "failed to allocate %zu input points", domain->num_points);
        return TM_STATUS_NO_MEMORY;
    }

    for (i = 0; i < domain->num_points; ++i) {
        tm_api_init_point(&points[i], domain->points[i].x, domain->points[i].y, (int) i, TM_VERTEX_INPUT);
    }

    *out_points = points;
    return TM_STATUS_OK;
}

static tm_status tm_copy_domain_to_pslg(const tm_domain *domain, TMPSLG *out_pslg, tm_error *out_error)
{
    tm_status public_status;
    size_t i;
    TMStatus status;

    memset(out_pslg, 0, sizeof(*out_pslg));

    public_status = tm_copy_points_to_internal(domain, &out_pslg->points, out_error);
    if (public_status != TM_STATUS_OK) {
        return public_status;
    }
    out_pslg->point_count = domain->num_points;

    if (domain->num_segments != 0) {
        out_pslg->segments = (TMSegment *) calloc(domain->num_segments, sizeof(*out_pslg->segments));
        if (out_pslg->segments == NULL) {
            tm_free_pslg(out_pslg);
            tm_api_set_error(out_error, TM_STATUS_NO_MEMORY, "failed to allocate %zu segments", domain->num_segments);
            return TM_STATUS_NO_MEMORY;
        }
        out_pslg->segment_count = domain->num_segments;

        for (i = 0; i < domain->num_segments; ++i) {
            if ((size_t) domain->segments[i].a >= domain->num_points ||
                (size_t) domain->segments[i].b >= domain->num_points) {
                tm_free_pslg(out_pslg);
                tm_api_set_error(
                    out_error,
                    TM_STATUS_GEOMETRY,
                    "segment %zu references point indices (%u, %u) outside 0..%zu",
                    i,
                    domain->segments[i].a,
                    domain->segments[i].b,
                    domain->num_points == 0 ? 0 : domain->num_points - 1
                );
                return TM_STATUS_GEOMETRY;
            }
            tm_api_init_segment(&out_pslg->segments[i], domain->segments[i].a, domain->segments[i].b, (int) i);
        }
    }

    if (domain->num_holes != 0) {
        out_pslg->holes = (double (*)[2]) calloc(domain->num_holes, sizeof(*out_pslg->holes));
        if (out_pslg->holes == NULL) {
            tm_free_pslg(out_pslg);
            tm_api_set_error(out_error, TM_STATUS_NO_MEMORY, "failed to allocate %zu hole markers", domain->num_holes);
            return TM_STATUS_NO_MEMORY;
        }
        out_pslg->hole_count = domain->num_holes;
        for (i = 0; i < domain->num_holes; ++i) {
            out_pslg->holes[i][0] = domain->holes[i].x;
            out_pslg->holes[i][1] = domain->holes[i].y;
        }
    }

    status = tm_validate_pslg(out_pslg);
    if (status != TM_OK) {
        tm_status mapped = tm_api_map_status(status);
        tm_free_pslg(out_pslg);
        tm_api_set_error(out_error, mapped, "%s", tm_internal_status_string(status));
        return mapped;
    }

    return TM_STATUS_OK;
}

static void tm_populate_public_mesh_metadata(const TMMesh *internal_mesh, tm_mesh *out_mesh)
{
    size_t point_index;
    size_t triangle_index;

    out_mesh->num_input_points = 0;
    out_mesh->num_segment_split_points = 0;
    out_mesh->num_triangle_split_points = 0;
    out_mesh->num_protected_corners = tm_count_protected_corners(internal_mesh);
    out_mesh->num_exempt_triangles = 0;

    for (point_index = 0; point_index < internal_mesh->point_count; ++point_index) {
        switch (internal_mesh->points[point_index].kind) {
            case TM_VERTEX_INPUT:
                out_mesh->num_input_points += 1;
                break;
            case TM_VERTEX_SEGMENT_SPLIT:
                out_mesh->num_segment_split_points += 1;
                break;
            case TM_VERTEX_TRIANGLE_SPLIT:
                out_mesh->num_triangle_split_points += 1;
                break;
            case TM_VERTEX_SUPER:
                break;
        }
    }

    for (triangle_index = 0; triangle_index < internal_mesh->triangle_count; ++triangle_index) {
        if (tm_triangle_is_quality_exempt(internal_mesh, (int) triangle_index)) {
            out_mesh->num_exempt_triangles += 1;
        }
    }
}

static tm_status tm_copy_internal_mesh(const TMMesh *internal_mesh, tm_mesh *out_mesh, tm_error *out_error)
{
    size_t i;
    size_t live_segment_count = 0;

    memset(out_mesh, 0, sizeof(*out_mesh));

    if (internal_mesh->point_count != 0) {
        out_mesh->points = (tm_point *) calloc(internal_mesh->point_count, sizeof(*out_mesh->points));
        if (out_mesh->points == NULL) {
            tm_api_set_error(out_error, TM_STATUS_NO_MEMORY, "failed to allocate %zu output points", internal_mesh->point_count);
            return TM_STATUS_NO_MEMORY;
        }
        out_mesh->num_points = internal_mesh->point_count;
        for (i = 0; i < internal_mesh->point_count; ++i) {
            out_mesh->points[i].x = internal_mesh->points[i].xy[0];
            out_mesh->points[i].y = internal_mesh->points[i].xy[1];
        }
    }

    if (internal_mesh->triangle_count != 0) {
        out_mesh->triangles = (uint32_t *) calloc(internal_mesh->triangle_count * 3, sizeof(*out_mesh->triangles));
        if (out_mesh->triangles == NULL) {
            tm_mesh_free(out_mesh);
            tm_api_set_error(
                out_error,
                TM_STATUS_NO_MEMORY,
                "failed to allocate %zu triangle indices",
                internal_mesh->triangle_count * 3
            );
            return TM_STATUS_NO_MEMORY;
        }
        out_mesh->num_triangles = internal_mesh->triangle_count;
        for (i = 0; i < internal_mesh->triangle_count; ++i) {
            out_mesh->triangles[i * 3 + 0] = (uint32_t) internal_mesh->triangles[i].v[0];
            out_mesh->triangles[i * 3 + 1] = (uint32_t) internal_mesh->triangles[i].v[1];
            out_mesh->triangles[i * 3 + 2] = (uint32_t) internal_mesh->triangles[i].v[2];
        }
    }

    for (i = 0; i < internal_mesh->segment_count; ++i) {
        if (internal_mesh->segments[i].live) {
            live_segment_count += 1;
        }
    }

    if (live_segment_count != 0) {
        size_t write_index = 0;

        out_mesh->segments = (tm_segment *) calloc(live_segment_count, sizeof(*out_mesh->segments));
        if (out_mesh->segments == NULL) {
            tm_mesh_free(out_mesh);
            tm_api_set_error(out_error, TM_STATUS_NO_MEMORY, "failed to allocate %zu output segments", live_segment_count);
            return TM_STATUS_NO_MEMORY;
        }
        out_mesh->num_segments = live_segment_count;

        for (i = 0; i < internal_mesh->segment_count; ++i) {
            if (!internal_mesh->segments[i].live) {
                continue;
            }
            out_mesh->segments[write_index].a = (uint32_t) internal_mesh->segments[i].v[0];
            out_mesh->segments[write_index].b = (uint32_t) internal_mesh->segments[i].v[1];
            write_index += 1;
        }
    }

    tm_populate_public_mesh_metadata(internal_mesh, out_mesh);
    return TM_STATUS_OK;
}

tm_status tm_api_map_status(TMStatus status)
{
    switch (status) {
        case TM_OK:
            return TM_STATUS_OK;
        case TM_ERR_ALLOC:
            return TM_STATUS_NO_MEMORY;
        case TM_ERR_IO:
            return TM_STATUS_IO;
        case TM_ERR_PARSE:
            return TM_STATUS_PARSE;
        case TM_ERR_DUPLICATE:
        case TM_ERR_TOO_FEW_POINTS:
        case TM_ERR_COLLINEAR:
        case TM_ERR_INVALID_PSLG:
        case TM_ERR_INVALID_MESH:
            return TM_STATUS_GEOMETRY;
        case TM_ERR_INTERNAL:
            return TM_STATUS_INTERNAL;
    }

    return TM_STATUS_INTERNAL;
}

void tm_api_clear_error(tm_error *error)
{
    if (error == NULL) {
        return;
    }

    error->code = TM_STATUS_OK;
    error->message[0] = '\0';
}

void tm_api_set_error_va(tm_error *error, tm_status status, const char *format, va_list args)
{
    if (error == NULL) {
        return;
    }

    error->code = status;
    if (format == NULL || format[0] == '\0') {
        error->message[0] = '\0';
        return;
    }

    vsnprintf(error->message, sizeof(error->message), format, args);
}

void tm_api_set_error(tm_error *error, tm_status status, const char *format, ...)
{
    va_list args;

    if (error == NULL) {
        return;
    }

    va_start(args, format);
    tm_api_set_error_va(error, status, format, args);
    va_end(args);
}

void tm_options_init(tm_options *options)
{
    if (options == NULL) {
        return;
    }

    options->min_angle_deg = 20.0;
    options->enable_refinement = 1;
    options->use_offcenters = 0;
    options->verbose = 0;
    options->enable_acute_protection = 1;
    options->acute_protection_mode = TM_ACUTE_PROTECTION_SHELL;
    options->protect_angle_deg = 0.0;
    options->max_refinement_steps = 0;
    options->max_protection_levels = 6;
}

const char *tm_status_string(tm_status status)
{
    switch (status) {
        case TM_STATUS_OK:
            return "ok";
        case TM_STATUS_INVALID_ARGUMENT:
            return "invalid argument";
        case TM_STATUS_GEOMETRY:
            return "invalid or unsupported geometry";
        case TM_STATUS_IO:
            return "I/O error";
        case TM_STATUS_PARSE:
            return "input parse error";
        case TM_STATUS_NO_MEMORY:
            return "allocation failed";
        case TM_STATUS_INTERNAL:
            return "internal error";
    }

    return "unknown status";
}

tm_status tm_generate(const tm_domain *domain, const tm_options *options, tm_mesh *out_mesh, tm_error *out_error)
{
    tm_options default_options;
    TMPoint *points = NULL;
    TMMesh internal_mesh;
    TMPSLG pslg;
    TMBuildOptions build_options;
    TMStatus internal_status;
    tm_status status;

    if (out_mesh == NULL) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "out_mesh must not be null");
        return TM_STATUS_INVALID_ARGUMENT;
    }

    memset(out_mesh, 0, sizeof(*out_mesh));
    memset(&internal_mesh, 0, sizeof(internal_mesh));
    memset(&pslg, 0, sizeof(pslg));
    tm_api_clear_error(out_error);

    status = tm_validate_domain_shape(domain, out_error);
    if (status != TM_STATUS_OK) {
        return status;
    }

    if (options == NULL) {
        tm_options_init(&default_options);
        options = &default_options;
    }

    if (options->min_angle_deg <= 0.0 || options->min_angle_deg >= 90.0) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "min_angle_deg must be in the range (0, 90)");
        return TM_STATUS_INVALID_ARGUMENT;
    }
    if (options->protect_angle_deg < 0.0 || options->protect_angle_deg >= 180.0) {
        tm_api_set_error(out_error, TM_STATUS_INVALID_ARGUMENT, "protect_angle_deg must be in the range [0, 180)");
        return TM_STATUS_INVALID_ARGUMENT;
    }

    build_options.verbose = options->verbose;
    build_options.refine = options->enable_refinement;
    build_options.use_offcenters = options->use_offcenters;
    build_options.protect_acute_corners = options->enable_acute_protection;
    build_options.acute_mode = options->acute_protection_mode == TM_ACUTE_PROTECTION_SIMPLE ?
        TM_ACUTE_MODE_SIMPLE :
        TM_ACUTE_MODE_SHELL;
    build_options.min_angle_deg = options->min_angle_deg;
    build_options.protect_angle_deg = options->protect_angle_deg;
    build_options.max_refinement_steps = options->max_refinement_steps;
    build_options.max_protection_levels = options->max_protection_levels;

    if (domain->num_segments == 0) {
        status = tm_copy_points_to_internal(domain, &points, out_error);
        if (status != TM_STATUS_OK) {
            return status;
        }

        internal_status = tm_build_mesh(points, domain->num_points, &internal_mesh);
        free(points);
        points = NULL;
    } else {
        status = tm_copy_domain_to_pslg(domain, &pslg, out_error);
        if (status != TM_STATUS_OK) {
            return status;
        }

        internal_status = tm_build_pslg_mesh(&pslg, &build_options, &internal_mesh);
        tm_free_pslg(&pslg);
    }

    if (internal_status != TM_OK) {
        status = tm_api_map_status(internal_status);
        tm_api_set_error(out_error, status, "%s", tm_internal_status_string(internal_status));
        return status;
    }

    status = tm_copy_internal_mesh(&internal_mesh, out_mesh, out_error);
    tm_free_mesh(&internal_mesh);
    return status;
}

void tm_mesh_free(tm_mesh *mesh)
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
    mesh->num_points = 0;
    mesh->num_triangles = 0;
    mesh->num_segments = 0;
    mesh->num_input_points = 0;
    mesh->num_segment_split_points = 0;
    mesh->num_triangle_split_points = 0;
    mesh->num_protected_corners = 0;
    mesh->num_exempt_triangles = 0;
}
