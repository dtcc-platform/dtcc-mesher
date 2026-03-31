#include <stdio.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"

int main(void)
{
    dtcc_mesher_point points[] = {
        {0.0, 0.0},
        {10.0, 0.0},
        {10.0, 10.0},
        {0.0, 10.0},
        {3.0, 3.0},
        {7.0, 3.0},
        {7.0, 7.0},
        {3.0, 7.0}
    };
    dtcc_mesher_segment segments[] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4}
    };
    dtcc_mesher_point holes[] = {
        {5.0, 5.0}
    };
    dtcc_mesher_domain domain;
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_error error;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_status status;

    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = segments;
    domain.num_segments = sizeof(segments) / sizeof(segments[0]);
    domain.holes = holes;
    domain.num_holes = sizeof(holes) / sizeof(holes[0]);

    dtcc_mesher_options_init(&options);
    options.verbose = 1;

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

    if (dtcc_mesher_write_svg(&mesh, "demo_basic.svg", &error) != DTCC_MESHER_STATUS_OK ||
        dtcc_mesher_write_quality_summary(&mesh, "demo_basic.summary.txt", &error) != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "write failed: %s\n", error.message[0] != '\0' ? error.message : "I/O error");
        dtcc_mesher_mesh_free(&mesh);
        return 1;
    }

    printf(
        "generated %zu points and %zu triangles (min angle %.2f deg)\n",
        mesh.num_points,
        mesh.num_triangles,
        summary.min_angle_deg_min
    );

    dtcc_mesher_mesh_free(&mesh);
    return 0;
}
