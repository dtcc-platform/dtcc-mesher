#ifndef DTCC_MESHER_VALIDATE_H
#define DTCC_MESHER_VALIDATE_H

#include "mesh.h"

typedef struct {
    size_t triangle_count;
    size_t boundary_edge_count;
    size_t adjacency_errors;
    size_t orientation_errors;
    size_t duplicate_triangle_errors;
    size_t incident_triangle_errors;
    size_t constrained_edge_errors;
    size_t local_delaunay_errors;
    size_t encroached_segment_errors;
    size_t bad_triangle_errors;
    size_t exempt_triangle_count;
    size_t first_encroached_segment;
    int first_encroaching_point;
    size_t first_bad_triangle;
} TMValidationReport;

TMStatus tm_validate_mesh(const TMMesh *mesh, int check_delaunay, TMValidationReport *out_report);
TMStatus tm_validate_quality_mesh(const TMMesh *mesh, double min_angle_deg, TMValidationReport *out_report);

#endif
