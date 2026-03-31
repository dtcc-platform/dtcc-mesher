#include <stdio.h>

#include "dtcc_mesher/dtcc_mesher.h"

int main(void)
{
    dtcc_mesher_point points[] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    dtcc_mesher_domain domain;
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_error error;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_status status;

    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = NULL;
    domain.num_segments = 0;
    domain.holes = NULL;
    domain.num_holes = 0;

    dtcc_mesher_options_init(&options);
    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "generate failed: %s\n", error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
        return 1;
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "analyze failed: %s\n", error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
        dtcc_mesher_mesh_free(&mesh);
        return 1;
    }

    printf("triangles: %zu\n", summary.triangle_count);
    printf("points: %zu\n", summary.point_count);
    printf("min angle: %.2f deg\n", summary.min_angle_deg_min);
    printf("max edge ratio: %.2f\n", summary.edge_ratio_max);

    dtcc_mesher_mesh_free(&mesh);
    return 0;
}
