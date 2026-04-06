#ifndef DTCC_MESHER_PUBLIC_DTCC_MESHER_H
#define DTCC_MESHER_PUBLIC_DTCC_MESHER_H

#include <stddef.h>
#include <stdint.h>

#include "dtcc_mesher_version.h"

#if defined(_WIN32) && defined(DTCC_MESHER_SHARED)
#  if defined(dtcc_mesher_EXPORTS)
#    define DTCC_MESHER_EXPORT __declspec(dllexport)
#  else
#    define DTCC_MESHER_EXPORT __declspec(dllimport)
#  endif
#else
#  define DTCC_MESHER_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x;
    double y;
} dtcc_mesher_point;

typedef struct {
    uint32_t a;
    uint32_t b;
} dtcc_mesher_segment;

typedef struct {
    dtcc_mesher_point *points;
    size_t num_points;
    dtcc_mesher_segment *segments;
    size_t num_segments;
    dtcc_mesher_point *holes;
    size_t num_holes;
} dtcc_mesher_domain;

typedef enum {
    DTCC_MESHER_ACUTE_PROTECTION_NONE = 0,
    DTCC_MESHER_ACUTE_PROTECTION_SIMPLE = 1,
    DTCC_MESHER_ACUTE_PROTECTION_SHELL = 2
} dtcc_mesher_acute_protection_mode;

typedef struct {
    double min_angle_deg;
    double max_area;
    double max_edge_length;
    int enable_refinement;
    int verbose;
    int enable_acute_protection;
    dtcc_mesher_acute_protection_mode acute_protection_mode;
    double protect_angle_deg;
    size_t max_refinement_steps;
    size_t max_protection_levels;
} dtcc_mesher_options;

typedef struct {
    dtcc_mesher_point *points;
    size_t num_points;

    uint32_t *triangles;
    size_t num_triangles;

    dtcc_mesher_segment *segments;
    size_t num_segments;

    size_t num_input_points;
    size_t num_segment_split_points;
    size_t num_triangle_split_points;
    size_t num_protected_corners;
    size_t num_exempt_triangles;
} dtcc_mesher_mesh;

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
} dtcc_mesher_quality_summary;

typedef enum {
    DTCC_MESHER_STATUS_OK = 0,
    DTCC_MESHER_STATUS_INVALID_ARGUMENT,
    DTCC_MESHER_STATUS_GEOMETRY,
    DTCC_MESHER_STATUS_IO,
    DTCC_MESHER_STATUS_PARSE,
    DTCC_MESHER_STATUS_NO_MEMORY,
    DTCC_MESHER_STATUS_INTERNAL
} dtcc_mesher_status;

typedef struct {
    dtcc_mesher_status code;
    char message[256];
} dtcc_mesher_error;

DTCC_MESHER_EXPORT void dtcc_mesher_options_init(dtcc_mesher_options *options);
DTCC_MESHER_EXPORT const char *dtcc_mesher_status_string(dtcc_mesher_status status);

DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_generate(
    const dtcc_mesher_domain *domain,
    const dtcc_mesher_options *options,
    dtcc_mesher_mesh *out_mesh,
    dtcc_mesher_error *out_error
);

DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_analyze_mesh(
    const dtcc_mesher_mesh *mesh,
    dtcc_mesher_quality_summary *out_summary,
    dtcc_mesher_error *out_error
);

DTCC_MESHER_EXPORT void dtcc_mesher_mesh_free(dtcc_mesher_mesh *mesh);

#ifdef __cplusplus
}
#endif

#endif
