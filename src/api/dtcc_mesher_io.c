#include "dtcc_mesher/dtcc_mesher_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_pslg.h"
#include "dtcc_mesher_api_internal.h"

static const double dtcc_mesher_api_pi = 3.14159265358979323846;

static const char *dtcc_mesher_internal_error_message(TMStatus status)
{
    const char *detail = NULL;

    if (status == TM_ERR_INVALID_PSLG) {
        detail = tm_last_pslg_error_detail();
        if (detail != NULL && detail[0] != '\0') {
            return detail;
        }
    }

    return tm_internal_status_string(status);
}

typedef struct {
    double area;
    double min_angle_deg;
    double max_angle_deg;
    double edge_ratio;
    double radius_edge_ratio;
} dtcc_mesher_triangle_metrics;

static int dtcc_mesher_has_suffix(const char *path, const char *suffix)
{
    size_t path_len;
    size_t suffix_len;

    if (path == NULL || suffix == NULL) {
        return 0;
    }

    path_len = strlen(path);
    suffix_len = strlen(suffix);
    return path_len >= suffix_len && strcmp(path + path_len - suffix_len, suffix) == 0;
}

static double dtcc_mesher_distance_xy(const dtcc_mesher_point *a, const dtcc_mesher_point *b)
{
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

static double dtcc_mesher_orient_xy(const dtcc_mesher_point *a, const dtcc_mesher_point *b, const dtcc_mesher_point *c)
{
    return (b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x);
}

static double dtcc_mesher_angle_from_lengths(double left, double right, double opposite)
{
    double cosine = (left * left + right * right - opposite * opposite) / (2.0 * left * right);

    if (cosine < -1.0) {
        cosine = -1.0;
    }
    if (cosine > 1.0) {
        cosine = 1.0;
    }

    return acos(cosine) * 180.0 / dtcc_mesher_api_pi;
}

static void dtcc_mesher_update_stats(double value, double *min_value, double *max_value, double *sum_value, size_t index)
{
    if (index == 0 || value < *min_value) {
        *min_value = value;
    }
    if (index == 0 || value > *max_value) {
        *max_value = value;
    }
    *sum_value += value;
}

static dtcc_mesher_status dtcc_mesher_validate_public_mesh(const dtcc_mesher_mesh *mesh, dtcc_mesher_error *out_error)
{
    size_t tri_index;

    if (mesh == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "mesh must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }
    if (mesh->num_points != 0 && mesh->points == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "mesh points must not be null when num_points > 0");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }
    if (mesh->num_triangles != 0 && mesh->triangles == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "mesh triangles must not be null when num_triangles > 0");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }
    if (mesh->num_segments != 0 && mesh->segments == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "mesh segments must not be null when num_segments > 0");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    for (tri_index = 0; tri_index < mesh->num_triangles; ++tri_index) {
        uint32_t a = mesh->triangles[tri_index * 3 + 0];
        uint32_t b = mesh->triangles[tri_index * 3 + 1];
        uint32_t c = mesh->triangles[tri_index * 3 + 2];

        if ((size_t) a >= mesh->num_points || (size_t) b >= mesh->num_points || (size_t) c >= mesh->num_points) {
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_GEOMETRY, "triangle %zu references an out-of-range point index", tri_index);
            return DTCC_MESHER_STATUS_GEOMETRY;
        }
    }

    return DTCC_MESHER_STATUS_OK;
}

static dtcc_mesher_status dtcc_mesher_compute_triangle_metrics(
    const dtcc_mesher_mesh *mesh,
    size_t tri_index,
    dtcc_mesher_triangle_metrics *out_metrics,
    dtcc_mesher_error *out_error
)
{
    const dtcc_mesher_point *a = &mesh->points[mesh->triangles[tri_index * 3 + 0]];
    const dtcc_mesher_point *b = &mesh->points[mesh->triangles[tri_index * 3 + 1]];
    const dtcc_mesher_point *c = &mesh->points[mesh->triangles[tri_index * 3 + 2]];
    double len_ab = dtcc_mesher_distance_xy(a, b);
    double len_bc = dtcc_mesher_distance_xy(b, c);
    double len_ca = dtcc_mesher_distance_xy(c, a);
    double shortest = len_ab;
    double longest = len_ab;
    double area = fabs(dtcc_mesher_orient_xy(a, b, c)) * 0.5;
    double angle_a;
    double angle_b;
    double angle_c;
    double min_angle;
    double max_angle;
    double edge_ratio;
    double circumradius;
    double radius_edge_ratio;

    if (len_bc < shortest) {
        shortest = len_bc;
    }
    if (len_ca < shortest) {
        shortest = len_ca;
    }
    if (len_bc > longest) {
        longest = len_bc;
    }
    if (len_ca > longest) {
        longest = len_ca;
    }

    if (area <= 0.0 || shortest <= 0.0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_GEOMETRY, "triangle %zu is degenerate", tri_index);
        return DTCC_MESHER_STATUS_GEOMETRY;
    }

    angle_a = dtcc_mesher_angle_from_lengths(len_ab, len_ca, len_bc);
    angle_b = dtcc_mesher_angle_from_lengths(len_ab, len_bc, len_ca);
    angle_c = dtcc_mesher_angle_from_lengths(len_bc, len_ca, len_ab);
    min_angle = angle_a;
    max_angle = angle_a;
    if (angle_b < min_angle) {
        min_angle = angle_b;
    }
    if (angle_c < min_angle) {
        min_angle = angle_c;
    }
    if (angle_b > max_angle) {
        max_angle = angle_b;
    }
    if (angle_c > max_angle) {
        max_angle = angle_c;
    }

    edge_ratio = longest / shortest;
    circumradius = (len_ab * len_bc * len_ca) / (4.0 * area);
    radius_edge_ratio = circumradius / shortest;

    out_metrics->area = area;
    out_metrics->min_angle_deg = min_angle;
    out_metrics->max_angle_deg = max_angle;
    out_metrics->edge_ratio = edge_ratio;
    out_metrics->radius_edge_ratio = radius_edge_ratio;
    return DTCC_MESHER_STATUS_OK;
}

dtcc_mesher_status dtcc_mesher_analyze_mesh(const dtcc_mesher_mesh *mesh, dtcc_mesher_quality_summary *out_summary, dtcc_mesher_error *out_error)
{
    double area_sum = 0.0;
    double min_angle_sum = 0.0;
    double edge_ratio_sum = 0.0;
    double radius_edge_ratio_sum = 0.0;
    size_t tri_index;
    dtcc_mesher_status status;

    dtcc_mesher_api_clear_error(out_error);
    if (out_summary == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "out_summary must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }
    memset(out_summary, 0, sizeof(*out_summary));

    status = dtcc_mesher_validate_public_mesh(mesh, out_error);
    if (status != DTCC_MESHER_STATUS_OK) {
        return status;
    }

    out_summary->point_count = mesh->num_points;
    out_summary->input_point_count = mesh->num_input_points;
    out_summary->segment_split_point_count = mesh->num_segment_split_points;
    out_summary->triangle_split_point_count = mesh->num_triangle_split_points;
    out_summary->steiner_point_count = mesh->num_segment_split_points + mesh->num_triangle_split_points;
    out_summary->protected_corner_count = mesh->num_protected_corners;
    out_summary->exempt_triangle_count = mesh->num_exempt_triangles;
    out_summary->triangle_count = mesh->num_triangles;

    for (tri_index = 0; tri_index < mesh->num_triangles; ++tri_index) {
        dtcc_mesher_triangle_metrics metrics;

        status = dtcc_mesher_compute_triangle_metrics(mesh, tri_index, &metrics, out_error);
        if (status != DTCC_MESHER_STATUS_OK) {
            return status;
        }

        dtcc_mesher_update_stats(metrics.area, &out_summary->area_min, &out_summary->area_max, &area_sum, tri_index);
        dtcc_mesher_update_stats(
            metrics.min_angle_deg,
            &out_summary->min_angle_deg_min,
            &out_summary->min_angle_deg_max,
            &min_angle_sum,
            tri_index
        );
        dtcc_mesher_update_stats(
            metrics.edge_ratio,
            &out_summary->edge_ratio_min,
            &out_summary->edge_ratio_max,
            &edge_ratio_sum,
            tri_index
        );
        dtcc_mesher_update_stats(
            metrics.radius_edge_ratio,
            &out_summary->radius_edge_ratio_min,
            &out_summary->radius_edge_ratio_max,
            &radius_edge_ratio_sum,
            tri_index
        );

        if (metrics.min_angle_deg < 20.0 - 1e-12) {
            out_summary->count_min_angle_lt_20 += 1;
        }
        if (metrics.min_angle_deg < 30.0 - 1e-12) {
            out_summary->count_min_angle_lt_30 += 1;
        }
    }

    if (mesh->num_triangles != 0) {
        out_summary->area_mean = area_sum / (double) mesh->num_triangles;
        out_summary->min_angle_deg_mean = min_angle_sum / (double) mesh->num_triangles;
        out_summary->edge_ratio_mean = edge_ratio_sum / (double) mesh->num_triangles;
        out_summary->radius_edge_ratio_mean = radius_edge_ratio_sum / (double) mesh->num_triangles;
    }

    return DTCC_MESHER_STATUS_OK;
}

dtcc_mesher_status dtcc_mesher_read_domain_file(const char *path, dtcc_mesher_domain *out_domain, dtcc_mesher_error *out_error)
{
    TMPoint *points = NULL;
    size_t point_count = 0;
    TMPSLG pslg;
    TMStatus internal_status;
    dtcc_mesher_status status;
    size_t i;

    dtcc_mesher_api_clear_error(out_error);
    if (path == NULL || out_domain == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "path and out_domain must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    memset(out_domain, 0, sizeof(*out_domain));
    memset(&pslg, 0, sizeof(pslg));

    if (dtcc_mesher_has_suffix(path, ".pslg")) {
        internal_status = tm_read_pslg_file(path, &pslg);
        if (internal_status != TM_OK) {
            status = dtcc_mesher_api_map_status(internal_status);
            dtcc_mesher_api_set_error(out_error, status, "%s", dtcc_mesher_internal_error_message(internal_status));
            return status;
        }

        out_domain->points = (dtcc_mesher_point *) calloc(pslg.point_count, sizeof(*out_domain->points));
        out_domain->segments = (dtcc_mesher_segment *) calloc(pslg.segment_count, sizeof(*out_domain->segments));
        out_domain->holes = (dtcc_mesher_point *) calloc(pslg.hole_count, sizeof(*out_domain->holes));
        if ((pslg.point_count != 0 && out_domain->points == NULL) ||
            (pslg.segment_count != 0 && out_domain->segments == NULL) ||
            (pslg.hole_count != 0 && out_domain->holes == NULL)) {
            tm_free_pslg(&pslg);
            dtcc_mesher_domain_free(out_domain);
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_NO_MEMORY, "failed to allocate domain arrays for %s", path);
            return DTCC_MESHER_STATUS_NO_MEMORY;
        }

        out_domain->num_points = pslg.point_count;
        out_domain->num_segments = pslg.segment_count;
        out_domain->num_holes = pslg.hole_count;

        for (i = 0; i < pslg.point_count; ++i) {
            out_domain->points[i].x = pslg.points[i].xy[0];
            out_domain->points[i].y = pslg.points[i].xy[1];
        }
        for (i = 0; i < pslg.segment_count; ++i) {
            out_domain->segments[i].a = (uint32_t) pslg.segments[i].v[0];
            out_domain->segments[i].b = (uint32_t) pslg.segments[i].v[1];
        }
        for (i = 0; i < pslg.hole_count; ++i) {
            out_domain->holes[i].x = pslg.holes[i][0];
            out_domain->holes[i].y = pslg.holes[i][1];
        }

        tm_free_pslg(&pslg);
        return DTCC_MESHER_STATUS_OK;
    }

    if (!dtcc_mesher_has_suffix(path, ".pts")) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "unsupported input extension for %s", path);
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    internal_status = tm_read_points_file(path, &points, &point_count);
    if (internal_status != TM_OK) {
        status = dtcc_mesher_api_map_status(internal_status);
        dtcc_mesher_api_set_error(out_error, status, "%s", dtcc_mesher_internal_error_message(internal_status));
        return status;
    }

    out_domain->points = (dtcc_mesher_point *) calloc(point_count, sizeof(*out_domain->points));
    if (out_domain->points == NULL && point_count != 0) {
        free(points);
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_NO_MEMORY, "failed to allocate %zu domain points", point_count);
        return DTCC_MESHER_STATUS_NO_MEMORY;
    }
    out_domain->num_points = point_count;

    for (i = 0; i < point_count; ++i) {
        out_domain->points[i].x = points[i].xy[0];
        out_domain->points[i].y = points[i].xy[1];
    }

    free(points);
    return DTCC_MESHER_STATUS_OK;
}

void dtcc_mesher_domain_free(dtcc_mesher_domain *domain)
{
    if (domain == NULL) {
        return;
    }

    free(domain->points);
    free(domain->segments);
    free(domain->holes);
    domain->points = NULL;
    domain->segments = NULL;
    domain->holes = NULL;
    domain->num_points = 0;
    domain->num_segments = 0;
    domain->num_holes = 0;
}

dtcc_mesher_status dtcc_mesher_write_triangles(const dtcc_mesher_mesh *mesh, const char *path, dtcc_mesher_error *out_error)
{
    FILE *stream;
    size_t tri_index;
    dtcc_mesher_status status;

    dtcc_mesher_api_clear_error(out_error);
    if (path == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "path must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    status = dtcc_mesher_validate_public_mesh(mesh, out_error);
    if (status != DTCC_MESHER_STATUS_OK) {
        return status;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to open %s for writing", path);
        return DTCC_MESHER_STATUS_IO;
    }

    for (tri_index = 0; tri_index < mesh->num_triangles; ++tri_index) {
        if (fprintf(
                stream,
                "%u %u %u\n",
                mesh->triangles[tri_index * 3 + 0],
                mesh->triangles[tri_index * 3 + 1],
                mesh->triangles[tri_index * 3 + 2]
            ) < 0) {
            fclose(stream);
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
            return DTCC_MESHER_STATUS_IO;
        }
    }

    if (fclose(stream) != 0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to close %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    return DTCC_MESHER_STATUS_OK;
}

dtcc_mesher_status dtcc_mesher_write_svg(const dtcc_mesher_mesh *mesh, const char *path, dtcc_mesher_error *out_error)
{
    FILE *stream;
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    double span_x;
    double span_y;
    double scale;
    double margin = 24.0;
    double canvas_w;
    double canvas_h;
    size_t i;
    dtcc_mesher_status status;

    dtcc_mesher_api_clear_error(out_error);
    if (path == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "path must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    status = dtcc_mesher_validate_public_mesh(mesh, out_error);
    if (status != DTCC_MESHER_STATUS_OK) {
        return status;
    }
    if (mesh->num_points == 0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "mesh must contain at least one point");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    min_x = mesh->points[0].x;
    min_y = mesh->points[0].y;
    max_x = mesh->points[0].x;
    max_y = mesh->points[0].y;

    for (i = 1; i < mesh->num_points; ++i) {
        if (mesh->points[i].x < min_x) {
            min_x = mesh->points[i].x;
        }
        if (mesh->points[i].x > max_x) {
            max_x = mesh->points[i].x;
        }
        if (mesh->points[i].y < min_y) {
            min_y = mesh->points[i].y;
        }
        if (mesh->points[i].y > max_y) {
            max_y = mesh->points[i].y;
        }
    }

    span_x = max_x - min_x;
    span_y = max_y - min_y;
    if (span_x == 0.0) {
        span_x = 1.0;
    }
    if (span_y == 0.0) {
        span_y = 1.0;
    }

    scale = 560.0 / ((span_x > span_y) ? span_x : span_y);
    canvas_w = span_x * scale + 2.0 * margin;
    canvas_h = span_y * scale + 2.0 * margin;

    stream = fopen(path, "w");
    if (stream == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to open %s for writing", path);
        return DTCC_MESHER_STATUS_IO;
    }

    if (fprintf(
            stream,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
            canvas_w,
            canvas_h,
            canvas_w,
            canvas_h
        ) < 0 ||
        fprintf(stream, "  <rect width=\"100%%\" height=\"100%%\" fill=\"white\"/>\n") < 0) {
        fclose(stream);
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    for (i = 0; i < mesh->num_triangles; ++i) {
        const dtcc_mesher_point *a = &mesh->points[mesh->triangles[i * 3 + 0]];
        const dtcc_mesher_point *b = &mesh->points[mesh->triangles[i * 3 + 1]];
        const dtcc_mesher_point *c = &mesh->points[mesh->triangles[i * 3 + 2]];
        double ax = margin + (a->x - min_x) * scale;
        double ay = canvas_h - margin - (a->y - min_y) * scale;
        double bx = margin + (b->x - min_x) * scale;
        double by = canvas_h - margin - (b->y - min_y) * scale;
        double cx = margin + (c->x - min_x) * scale;
        double cy = canvas_h - margin - (c->y - min_y) * scale;

        if (fprintf(
                stream,
                "  <polygon points=\"%.6f,%.6f %.6f,%.6f %.6f,%.6f\" fill=\"none\" stroke=\"black\" stroke-width=\"1.25\"/>\n",
                ax,
                ay,
                bx,
                by,
                cx,
                cy
            ) < 0) {
            fclose(stream);
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
            return DTCC_MESHER_STATUS_IO;
        }
    }

    for (i = 0; i < mesh->num_segments; ++i) {
        const dtcc_mesher_point *a = &mesh->points[mesh->segments[i].a];
        const dtcc_mesher_point *b = &mesh->points[mesh->segments[i].b];
        double ax = margin + (a->x - min_x) * scale;
        double ay = canvas_h - margin - (a->y - min_y) * scale;
        double bx = margin + (b->x - min_x) * scale;
        double by = canvas_h - margin - (b->y - min_y) * scale;

        if (fprintf(
                stream,
                "  <line x1=\"%.6f\" y1=\"%.6f\" x2=\"%.6f\" y2=\"%.6f\" stroke=\"#b22222\" stroke-width=\"3\"/>\n",
                ax,
                ay,
                bx,
                by
            ) < 0) {
            fclose(stream);
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
            return DTCC_MESHER_STATUS_IO;
        }
    }

    if (mesh->num_points <= 16) {
        for (i = 0; i < mesh->num_points; ++i) {
            double px = margin + (mesh->points[i].x - min_x) * scale;
            double py = canvas_h - margin - (mesh->points[i].y - min_y) * scale;

            if (fprintf(
                    stream,
                    "  <text x=\"%.6f\" y=\"%.6f\" font-family=\"monospace\" font-size=\"12\" fill=\"black\">%zu</text>\n",
                    px + 5.0,
                    py - 5.0,
                    i
                ) < 0) {
                fclose(stream);
                dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
                return DTCC_MESHER_STATUS_IO;
            }
        }
    }

    if (fprintf(stream, "</svg>\n") < 0 || fclose(stream) != 0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to finalize %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    return DTCC_MESHER_STATUS_OK;
}

dtcc_mesher_status dtcc_mesher_write_quality_csv(const dtcc_mesher_mesh *mesh, const char *path, dtcc_mesher_error *out_error)
{
    FILE *stream;
    size_t tri_index;
    dtcc_mesher_status status;

    dtcc_mesher_api_clear_error(out_error);
    if (path == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "path must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    status = dtcc_mesher_validate_public_mesh(mesh, out_error);
    if (status != DTCC_MESHER_STATUS_OK) {
        return status;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to open %s for writing", path);
        return DTCC_MESHER_STATUS_IO;
    }

    if (fprintf(
            stream,
            "tri_id,v0,v1,v2,area,min_angle_deg,max_angle_deg,edge_ratio,radius_edge_ratio\n"
        ) < 0) {
        fclose(stream);
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    for (tri_index = 0; tri_index < mesh->num_triangles; ++tri_index) {
        dtcc_mesher_triangle_metrics metrics;

        status = dtcc_mesher_compute_triangle_metrics(mesh, tri_index, &metrics, out_error);
        if (status != DTCC_MESHER_STATUS_OK) {
            fclose(stream);
            return status;
        }

        if (fprintf(
                stream,
                "%zu,%u,%u,%u,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                tri_index,
                mesh->triangles[tri_index * 3 + 0],
                mesh->triangles[tri_index * 3 + 1],
                mesh->triangles[tri_index * 3 + 2],
                metrics.area,
                metrics.min_angle_deg,
                metrics.max_angle_deg,
                metrics.edge_ratio,
                metrics.radius_edge_ratio
            ) < 0) {
            fclose(stream);
            dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
            return DTCC_MESHER_STATUS_IO;
        }
    }

    if (fclose(stream) != 0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to close %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    return DTCC_MESHER_STATUS_OK;
}

dtcc_mesher_status dtcc_mesher_write_quality_summary(const dtcc_mesher_mesh *mesh, const char *path, dtcc_mesher_error *out_error)
{
    FILE *stream;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_status status;

    dtcc_mesher_api_clear_error(out_error);
    if (path == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_INVALID_ARGUMENT, "path must not be null");
        return DTCC_MESHER_STATUS_INVALID_ARGUMENT;
    }

    status = dtcc_mesher_analyze_mesh(mesh, &summary, out_error);
    if (status != DTCC_MESHER_STATUS_OK) {
        return status;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to open %s for writing", path);
        return DTCC_MESHER_STATUS_IO;
    }

    if (fprintf(stream, "triangle_count=%zu\n", summary.triangle_count) < 0 ||
        fprintf(stream, "point_count=%zu\n", summary.point_count) < 0 ||
        fprintf(stream, "input_point_count=%zu\n", summary.input_point_count) < 0 ||
        fprintf(stream, "steiner_point_count=%zu\n", summary.steiner_point_count) < 0 ||
        fprintf(stream, "segment_split_point_count=%zu\n", summary.segment_split_point_count) < 0 ||
        fprintf(stream, "triangle_split_point_count=%zu\n", summary.triangle_split_point_count) < 0 ||
        fprintf(stream, "protected_corner_count=%zu\n", summary.protected_corner_count) < 0 ||
        fprintf(stream, "exempt_triangle_count=%zu\n", summary.exempt_triangle_count) < 0 ||
        fprintf(stream, "area_min=%.17g\n", summary.area_min) < 0 ||
        fprintf(stream, "area_mean=%.17g\n", summary.area_mean) < 0 ||
        fprintf(stream, "area_max=%.17g\n", summary.area_max) < 0 ||
        fprintf(stream, "min_angle_deg_min=%.17g\n", summary.min_angle_deg_min) < 0 ||
        fprintf(stream, "min_angle_deg_mean=%.17g\n", summary.min_angle_deg_mean) < 0 ||
        fprintf(stream, "min_angle_deg_max=%.17g\n", summary.min_angle_deg_max) < 0 ||
        fprintf(stream, "edge_ratio_min=%.17g\n", summary.edge_ratio_min) < 0 ||
        fprintf(stream, "edge_ratio_mean=%.17g\n", summary.edge_ratio_mean) < 0 ||
        fprintf(stream, "edge_ratio_max=%.17g\n", summary.edge_ratio_max) < 0 ||
        fprintf(stream, "radius_edge_ratio_min=%.17g\n", summary.radius_edge_ratio_min) < 0 ||
        fprintf(stream, "radius_edge_ratio_mean=%.17g\n", summary.radius_edge_ratio_mean) < 0 ||
        fprintf(stream, "radius_edge_ratio_max=%.17g\n", summary.radius_edge_ratio_max) < 0 ||
        fprintf(stream, "count_min_angle_lt_20=%zu\n", summary.count_min_angle_lt_20) < 0 ||
        fprintf(stream, "count_min_angle_lt_30=%zu\n", summary.count_min_angle_lt_30) < 0) {
        fclose(stream);
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to write %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    if (fclose(stream) != 0) {
        dtcc_mesher_api_set_error(out_error, DTCC_MESHER_STATUS_IO, "failed to close %s", path);
        return DTCC_MESHER_STATUS_IO;
    }

    return DTCC_MESHER_STATUS_OK;
}
