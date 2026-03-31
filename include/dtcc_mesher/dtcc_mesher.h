#ifndef DTCC_MESHER_PUBLIC_DTCC_MESHER_H
#define DTCC_MESHER_PUBLIC_DTCC_MESHER_H

#include <stddef.h>
#include <stdint.h>

#include "dtcc_mesher_version.h"

#if defined(_WIN32) && defined(DTCC_MESHER_SHARED)
#  if defined(dtcc_mesher_EXPORTS)
#    define TM_EXPORT __declspec(dllexport)
#  else
#    define TM_EXPORT __declspec(dllimport)
#  endif
#else
#  define TM_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x;
    double y;
} tm_point;

typedef struct {
    uint32_t a;
    uint32_t b;
} tm_segment;

typedef struct {
    tm_point *points;
    size_t num_points;
    tm_segment *segments;
    size_t num_segments;
    tm_point *holes;
    size_t num_holes;
} tm_domain;

typedef enum {
    TM_ACUTE_PROTECTION_NONE = 0,
    TM_ACUTE_PROTECTION_SIMPLE = 1,
    TM_ACUTE_PROTECTION_SHELL = 2
} tm_acute_protection_mode;

typedef struct {
    double min_angle_deg;
    int enable_refinement;
    int use_offcenters;
    int verbose;
    int enable_acute_protection;
    tm_acute_protection_mode acute_protection_mode;
    double protect_angle_deg;
    size_t max_refinement_steps;
    size_t max_protection_levels;
} tm_options;

typedef struct {
    tm_point *points;
    size_t num_points;

    uint32_t *triangles;
    size_t num_triangles;

    tm_segment *segments;
    size_t num_segments;

    size_t num_input_points;
    size_t num_segment_split_points;
    size_t num_triangle_split_points;
    size_t num_protected_corners;
    size_t num_exempt_triangles;
} tm_mesh;

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
} tm_quality_summary;

typedef enum {
    TM_STATUS_OK = 0,
    TM_STATUS_INVALID_ARGUMENT,
    TM_STATUS_GEOMETRY,
    TM_STATUS_IO,
    TM_STATUS_PARSE,
    TM_STATUS_NO_MEMORY,
    TM_STATUS_INTERNAL
} tm_status;

typedef struct {
    tm_status code;
    char message[256];
} tm_error;

TM_EXPORT void tm_options_init(tm_options *options);
TM_EXPORT const char *tm_status_string(tm_status status);

TM_EXPORT tm_status tm_generate(
    const tm_domain *domain,
    const tm_options *options,
    tm_mesh *out_mesh,
    tm_error *out_error
);

TM_EXPORT tm_status tm_analyze_mesh(
    const tm_mesh *mesh,
    tm_quality_summary *out_summary,
    tm_error *out_error
);

TM_EXPORT void tm_mesh_free(tm_mesh *mesh);

#ifdef __cplusplus
}
#endif

#endif
