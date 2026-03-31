#include "report.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "predicates.h"

static const double tm_pi = 3.14159265358979323846;

static double tm_distance(const TMPoint *a, const TMPoint *b)
{
    double dx = a->xy[0] - b->xy[0];
    double dy = a->xy[1] - b->xy[1];
    return sqrt(dx * dx + dy * dy);
}

static double tm_angle_from_lengths(double left, double right, double opposite)
{
    double cosine = (left * left + right * right - opposite * opposite) / (2.0 * left * right);

    if (cosine < -1.0) {
        cosine = -1.0;
    }
    if (cosine > 1.0) {
        cosine = 1.0;
    }

    return acos(cosine) * 180.0 / tm_pi;
}

static void tm_update_stats(double value, double *min_value, double *max_value, double *sum_value, size_t index)
{
    if (index == 0 || value < *min_value) {
        *min_value = value;
    }
    if (index == 0 || value > *max_value) {
        *max_value = value;
    }
    *sum_value += value;
}

TMStatus tm_compute_metrics(const TMMesh *mesh, TMTriangleMetric **out_metrics, TMSummary *out_summary)
{
    TMTriangleMetric *metrics = NULL;
    double area_sum = 0.0;
    double min_angle_sum = 0.0;
    double edge_ratio_sum = 0.0;
    double radius_edge_ratio_sum = 0.0;
    size_t i;

    if (mesh == NULL || out_metrics == NULL || out_summary == NULL) {
        return TM_ERR_INTERNAL;
    }

    *out_metrics = NULL;
    out_summary->point_count = mesh->point_count;
    out_summary->input_point_count = 0;
    out_summary->steiner_point_count = 0;
    out_summary->segment_split_point_count = 0;
    out_summary->triangle_split_point_count = 0;
    out_summary->protected_corner_count = tm_count_protected_corners(mesh);
    out_summary->exempt_triangle_count = 0;
    out_summary->triangle_count = mesh->triangle_count;
    out_summary->area_min = 0.0;
    out_summary->area_mean = 0.0;
    out_summary->area_max = 0.0;
    out_summary->min_angle_deg_min = 0.0;
    out_summary->min_angle_deg_mean = 0.0;
    out_summary->min_angle_deg_max = 0.0;
    out_summary->edge_ratio_min = 0.0;
    out_summary->edge_ratio_mean = 0.0;
    out_summary->edge_ratio_max = 0.0;
    out_summary->radius_edge_ratio_min = 0.0;
    out_summary->radius_edge_ratio_mean = 0.0;
    out_summary->radius_edge_ratio_max = 0.0;
    out_summary->count_min_angle_lt_20 = 0;
    out_summary->count_min_angle_lt_30 = 0;

    for (i = 0; i < mesh->point_count; ++i) {
        switch (mesh->points[i].kind) {
            case TM_VERTEX_INPUT:
                out_summary->input_point_count += 1;
                break;
            case TM_VERTEX_SEGMENT_SPLIT:
                out_summary->segment_split_point_count += 1;
                out_summary->steiner_point_count += 1;
                break;
            case TM_VERTEX_TRIANGLE_SPLIT:
                out_summary->triangle_split_point_count += 1;
                out_summary->steiner_point_count += 1;
                break;
            case TM_VERTEX_SUPER:
                break;
        }
    }

    if (mesh->triangle_count == 0) {
        return TM_OK;
    }

    metrics = (TMTriangleMetric *) calloc(mesh->triangle_count, sizeof(*metrics));
    if (metrics == NULL) {
        return TM_ERR_ALLOC;
    }

    for (i = 0; i < mesh->triangle_count; ++i) {
        const TMTriangle *triangle = &mesh->triangles[i];
        const TMPoint *a = &mesh->points[triangle->v[0]];
        const TMPoint *b = &mesh->points[triangle->v[1]];
        const TMPoint *c = &mesh->points[triangle->v[2]];
        double len_ab = tm_distance(a, b);
        double len_bc = tm_distance(b, c);
        double len_ca = tm_distance(c, a);
        double shortest = len_ab;
        double longest = len_ab;
        double area = fabs(orient2d((REAL *) a->xy, (REAL *) b->xy, (REAL *) c->xy)) * 0.5;
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
            free(metrics);
            return TM_ERR_INTERNAL;
        }

        angle_a = tm_angle_from_lengths(len_ab, len_ca, len_bc);
        angle_b = tm_angle_from_lengths(len_ab, len_bc, len_ca);
        angle_c = tm_angle_from_lengths(len_bc, len_ca, len_ab);
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

        metrics[i].tri_id = i;
        metrics[i].v0 = triangle->v[0];
        metrics[i].v1 = triangle->v[1];
        metrics[i].v2 = triangle->v[2];
        metrics[i].area = area;
        metrics[i].min_angle_deg = min_angle;
        metrics[i].max_angle_deg = max_angle;
        metrics[i].edge_ratio = edge_ratio;
        metrics[i].radius_edge_ratio = radius_edge_ratio;

        if (tm_triangle_is_quality_exempt(mesh, (int) i)) {
            out_summary->exempt_triangle_count += 1;
        }

        tm_update_stats(area, &out_summary->area_min, &out_summary->area_max, &area_sum, i);
        tm_update_stats(min_angle, &out_summary->min_angle_deg_min, &out_summary->min_angle_deg_max, &min_angle_sum, i);
        tm_update_stats(edge_ratio, &out_summary->edge_ratio_min, &out_summary->edge_ratio_max, &edge_ratio_sum, i);
        tm_update_stats(
            radius_edge_ratio,
            &out_summary->radius_edge_ratio_min,
            &out_summary->radius_edge_ratio_max,
            &radius_edge_ratio_sum,
            i
        );

        if (min_angle < 20.0 - 1e-12) {
            out_summary->count_min_angle_lt_20 += 1;
        }
        if (min_angle < 30.0 - 1e-12) {
            out_summary->count_min_angle_lt_30 += 1;
        }
    }

    out_summary->area_mean = area_sum / (double) mesh->triangle_count;
    out_summary->min_angle_deg_mean = min_angle_sum / (double) mesh->triangle_count;
    out_summary->edge_ratio_mean = edge_ratio_sum / (double) mesh->triangle_count;
    out_summary->radius_edge_ratio_mean = radius_edge_ratio_sum / (double) mesh->triangle_count;
    *out_metrics = metrics;
    return TM_OK;
}

void tm_free_metrics(TMTriangleMetric *metrics)
{
    free(metrics);
}

TMStatus tm_write_tri_file(const char *path, const TMMesh *mesh)
{
    FILE *stream;
    size_t i;

    if (path == NULL || mesh == NULL) {
        return TM_ERR_INTERNAL;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        return TM_ERR_IO;
    }

    for (i = 0; i < mesh->triangle_count; ++i) {
        if (fprintf(
                stream,
                "%d %d %d\n",
                mesh->triangles[i].v[0],
                mesh->triangles[i].v[1],
                mesh->triangles[i].v[2]
            ) < 0) {
            fclose(stream);
            return TM_ERR_IO;
        }
    }

    if (fclose(stream) != 0) {
        return TM_ERR_IO;
    }

    return TM_OK;
}

TMStatus tm_write_svg_file(const char *path, const TMMesh *mesh)
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

    if (path == NULL || mesh == NULL || mesh->point_count == 0) {
        return TM_ERR_INTERNAL;
    }

    min_x = mesh->points[0].xy[0];
    min_y = mesh->points[0].xy[1];
    max_x = mesh->points[0].xy[0];
    max_y = mesh->points[0].xy[1];

    for (i = 1; i < mesh->point_count; ++i) {
        if (mesh->points[i].xy[0] < min_x) {
            min_x = mesh->points[i].xy[0];
        }
        if (mesh->points[i].xy[0] > max_x) {
            max_x = mesh->points[i].xy[0];
        }
        if (mesh->points[i].xy[1] < min_y) {
            min_y = mesh->points[i].xy[1];
        }
        if (mesh->points[i].xy[1] > max_y) {
            max_y = mesh->points[i].xy[1];
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
        return TM_ERR_IO;
    }

    if (fprintf(
            stream,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
            canvas_w,
            canvas_h,
            canvas_w,
            canvas_h
        ) < 0) {
        fclose(stream);
        return TM_ERR_IO;
    }

    if (fprintf(stream, "  <rect width=\"100%%\" height=\"100%%\" fill=\"white\"/>\n") < 0) {
        fclose(stream);
        return TM_ERR_IO;
    }

    for (i = 0; i < mesh->triangle_count; ++i) {
        const TMTriangle *triangle = &mesh->triangles[i];
        const TMPoint *a = &mesh->points[triangle->v[0]];
        const TMPoint *b = &mesh->points[triangle->v[1]];
        const TMPoint *c = &mesh->points[triangle->v[2]];
        double ax = margin + (a->xy[0] - min_x) * scale;
        double ay = canvas_h - margin - (a->xy[1] - min_y) * scale;
        double bx = margin + (b->xy[0] - min_x) * scale;
        double by = canvas_h - margin - (b->xy[1] - min_y) * scale;
        double cx = margin + (c->xy[0] - min_x) * scale;
        double cy = canvas_h - margin - (c->xy[1] - min_y) * scale;

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
            return TM_ERR_IO;
        }
    }

    for (i = 0; i < mesh->triangle_count; ++i) {
        const TMTriangle *triangle = &mesh->triangles[i];
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            int a_index;
            int b_index;
            int neighbor = triangle->nbr[edge];
            const TMPoint *a;
            const TMPoint *b;
            double ax;
            double ay;
            double bx;
            double by;

            if (!triangle->constrained[edge]) {
                continue;
            }
            if (neighbor >= 0 && (size_t) neighbor < i) {
                continue;
            }

            tm_triangle_edge_vertices(triangle, edge, &a_index, &b_index);
            a = &mesh->points[a_index];
            b = &mesh->points[b_index];
            ax = margin + (a->xy[0] - min_x) * scale;
            ay = canvas_h - margin - (a->xy[1] - min_y) * scale;
            bx = margin + (b->xy[0] - min_x) * scale;
            by = canvas_h - margin - (b->xy[1] - min_y) * scale;

            if (fprintf(
                    stream,
                    "  <line x1=\"%.6f\" y1=\"%.6f\" x2=\"%.6f\" y2=\"%.6f\" stroke=\"#b22222\" stroke-width=\"3\"/>\n",
                    ax,
                    ay,
                    bx,
                    by
                ) < 0) {
                fclose(stream);
                return TM_ERR_IO;
            }
        }
    }

    for (i = 0; i < mesh->point_count; ++i) {
        double px = margin + (mesh->points[i].xy[0] - min_x) * scale;
        double py = canvas_h - margin - (mesh->points[i].xy[1] - min_y) * scale;

        if (mesh->point_count <= 16) {
            if (fprintf(
                    stream,
                    "  <text x=\"%.6f\" y=\"%.6f\" font-family=\"monospace\" font-size=\"12\" fill=\"black\">%zu</text>\n",
                    px + 5.0,
                    py - 5.0,
                    i
                ) < 0) {
                fclose(stream);
                return TM_ERR_IO;
            }
        }
    }

    if (fprintf(stream, "</svg>\n") < 0) {
        fclose(stream);
        return TM_ERR_IO;
    }

    if (fclose(stream) != 0) {
        return TM_ERR_IO;
    }

    return TM_OK;
}

TMStatus tm_write_metrics_csv(const char *path, const TMTriangleMetric *metrics, size_t metric_count)
{
    FILE *stream;
    size_t i;

    if (path == NULL || (metrics == NULL && metric_count != 0)) {
        return TM_ERR_INTERNAL;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        return TM_ERR_IO;
    }

    if (fprintf(
            stream,
            "tri_id,v0,v1,v2,area,min_angle_deg,max_angle_deg,edge_ratio,radius_edge_ratio\n"
        ) < 0) {
        fclose(stream);
        return TM_ERR_IO;
    }

    for (i = 0; i < metric_count; ++i) {
        if (fprintf(
                stream,
                "%zu,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                metrics[i].tri_id,
                metrics[i].v0,
                metrics[i].v1,
                metrics[i].v2,
                metrics[i].area,
                metrics[i].min_angle_deg,
                metrics[i].max_angle_deg,
                metrics[i].edge_ratio,
                metrics[i].radius_edge_ratio
            ) < 0) {
            fclose(stream);
            return TM_ERR_IO;
        }
    }

    if (fclose(stream) != 0) {
        return TM_ERR_IO;
    }

    return TM_OK;
}

TMStatus tm_write_summary_file(const char *path, const TMSummary *summary)
{
    FILE *stream;

    if (path == NULL || summary == NULL) {
        return TM_ERR_INTERNAL;
    }

    stream = fopen(path, "w");
    if (stream == NULL) {
        return TM_ERR_IO;
    }

    if (fprintf(stream, "triangle_count=%zu\n", summary->triangle_count) < 0 ||
        fprintf(stream, "point_count=%zu\n", summary->point_count) < 0 ||
        fprintf(stream, "input_point_count=%zu\n", summary->input_point_count) < 0 ||
        fprintf(stream, "steiner_point_count=%zu\n", summary->steiner_point_count) < 0 ||
        fprintf(stream, "segment_split_point_count=%zu\n", summary->segment_split_point_count) < 0 ||
        fprintf(stream, "triangle_split_point_count=%zu\n", summary->triangle_split_point_count) < 0 ||
        fprintf(stream, "protected_corner_count=%zu\n", summary->protected_corner_count) < 0 ||
        fprintf(stream, "exempt_triangle_count=%zu\n", summary->exempt_triangle_count) < 0 ||
        fprintf(stream, "area_min=%.17g\n", summary->area_min) < 0 ||
        fprintf(stream, "area_mean=%.17g\n", summary->area_mean) < 0 ||
        fprintf(stream, "area_max=%.17g\n", summary->area_max) < 0 ||
        fprintf(stream, "min_angle_deg_min=%.17g\n", summary->min_angle_deg_min) < 0 ||
        fprintf(stream, "min_angle_deg_mean=%.17g\n", summary->min_angle_deg_mean) < 0 ||
        fprintf(stream, "min_angle_deg_max=%.17g\n", summary->min_angle_deg_max) < 0 ||
        fprintf(stream, "edge_ratio_min=%.17g\n", summary->edge_ratio_min) < 0 ||
        fprintf(stream, "edge_ratio_mean=%.17g\n", summary->edge_ratio_mean) < 0 ||
        fprintf(stream, "edge_ratio_max=%.17g\n", summary->edge_ratio_max) < 0 ||
        fprintf(stream, "radius_edge_ratio_min=%.17g\n", summary->radius_edge_ratio_min) < 0 ||
        fprintf(stream, "radius_edge_ratio_mean=%.17g\n", summary->radius_edge_ratio_mean) < 0 ||
        fprintf(stream, "radius_edge_ratio_max=%.17g\n", summary->radius_edge_ratio_max) < 0 ||
        fprintf(stream, "count_min_angle_lt_20=%zu\n", summary->count_min_angle_lt_20) < 0 ||
        fprintf(stream, "count_min_angle_lt_30=%zu\n", summary->count_min_angle_lt_30) < 0) {
        fclose(stream);
        return TM_ERR_IO;
    }

    if (fclose(stream) != 0) {
        return TM_ERR_IO;
    }

    return TM_OK;
}
