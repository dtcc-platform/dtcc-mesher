#include <stdio.h>

#include "dtcc_mesher/dtcc_mesher.h"

int main(void)
{
    tm_point points[] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    tm_domain domain;
    tm_options options;
    tm_mesh mesh;
    tm_error error;
    tm_quality_summary summary;
    tm_status status;

    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = NULL;
    domain.num_segments = 0;
    domain.holes = NULL;
    domain.num_holes = 0;

    tm_options_init(&options);
    status = tm_generate(&domain, &options, &mesh, &error);
    if (status != TM_STATUS_OK) {
        fprintf(stderr, "generate failed: %s\n", error.message[0] != '\0' ? error.message : tm_status_string(status));
        return 1;
    }

    status = tm_analyze_mesh(&mesh, &summary, &error);
    if (status != TM_STATUS_OK) {
        fprintf(stderr, "analyze failed: %s\n", error.message[0] != '\0' ? error.message : tm_status_string(status));
        tm_mesh_free(&mesh);
        return 1;
    }

    printf("triangles: %zu\n", summary.triangle_count);
    printf("points: %zu\n", summary.point_count);
    printf("min angle: %.2f deg\n", summary.min_angle_deg_min);
    printf("max edge ratio: %.2f\n", summary.edge_ratio_max);

    tm_mesh_free(&mesh);
    return 0;
}
