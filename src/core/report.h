#ifndef DTCC_MESHER_REPORT_H
#define DTCC_MESHER_REPORT_H

#include <stddef.h>

#include "mesh.h"

typedef struct {
    size_t tri_id;
    int v0;
    int v1;
    int v2;
    double area;
    double min_angle_deg;
    double max_angle_deg;
    double edge_ratio;
    double radius_edge_ratio;
} TMTriangleMetric;

typedef struct {
    size_t point_count;
    size_t input_point_count;
    size_t steiner_point_count;
    size_t segment_split_point_count;
    size_t triangle_split_point_count;
    size_t protected_corner_count;
    size_t exempt_triangle_count;
    size_t triangle_count;
    double area_min;
    double area_mean;
    double area_max;
    double min_angle_deg_min;
    double min_angle_deg_mean;
    double min_angle_deg_max;
    double edge_ratio_min;
    double edge_ratio_mean;
    double edge_ratio_max;
    double radius_edge_ratio_min;
    double radius_edge_ratio_mean;
    double radius_edge_ratio_max;
    size_t count_min_angle_lt_20;
    size_t count_min_angle_lt_30;
} TMSummary;

TMStatus tm_compute_metrics(const TMMesh *mesh, TMTriangleMetric **out_metrics, TMSummary *out_summary);
void tm_free_metrics(TMTriangleMetric *metrics);
TMStatus tm_write_tri_file(const char *path, const TMMesh *mesh);
TMStatus tm_write_svg_file(const char *path, const TMMesh *mesh);
TMStatus tm_write_metrics_csv(const char *path, const TMTriangleMetric *metrics, size_t metric_count);
TMStatus tm_write_summary_file(const char *path, const TMSummary *summary);

#endif
