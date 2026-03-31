#include <stdio.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"

int main(void)
{
    tm_point points[] = {
        {0.0, 0.0},
        {10.0, 0.0},
        {10.0, 10.0},
        {0.0, 10.0},
        {3.0, 3.0},
        {7.0, 3.0},
        {7.0, 7.0},
        {3.0, 7.0}
    };
    tm_segment segments[] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4}
    };
    tm_point holes[] = {
        {5.0, 5.0}
    };
    tm_domain domain;
    tm_options options;
    tm_mesh mesh;
    tm_error error;
    tm_quality_summary summary;
    tm_status status;

    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = segments;
    domain.num_segments = sizeof(segments) / sizeof(segments[0]);
    domain.holes = holes;
    domain.num_holes = sizeof(holes) / sizeof(holes[0]);

    tm_options_init(&options);
    options.verbose = 1;

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

    if (tm_write_svg(&mesh, "demo_basic.svg", &error) != TM_STATUS_OK ||
        tm_write_quality_summary(&mesh, "demo_basic.summary.txt", &error) != TM_STATUS_OK) {
        fprintf(stderr, "write failed: %s\n", error.message[0] != '\0' ? error.message : "I/O error");
        tm_mesh_free(&mesh);
        return 1;
    }

    printf(
        "generated %zu points and %zu triangles (min angle %.2f deg)\n",
        mesh.num_points,
        mesh.num_triangles,
        summary.min_angle_deg_min
    );

    tm_mesh_free(&mesh);
    return 0;
}
