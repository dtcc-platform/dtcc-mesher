#include <stdio.h>
#include <string.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"

int main(int argc, char **argv)
{
    tm_domain domain;
    tm_options options;
    tm_mesh mesh;
    tm_quality_summary summary;
    tm_error error;
    tm_status status;

    if (argc != 2) {
        fprintf(stderr, "usage: api_dump input.(pts|pslg)\n");
        return 2;
    }

    memset(&domain, 0, sizeof(domain));
    memset(&mesh, 0, sizeof(mesh));
    memset(&summary, 0, sizeof(summary));
    memset(&error, 0, sizeof(error));

    status = tm_read_domain_file(argv[1], &domain, &error);
    if (status != TM_STATUS_OK) {
        fprintf(stderr, "read failed: %s\n", error.message[0] != '\0' ? error.message : tm_status_string(status));
        return 1;
    }

    tm_options_init(&options);
    status = tm_generate(&domain, &options, &mesh, &error);
    if (status != TM_STATUS_OK) {
        fprintf(stderr, "generate failed: %s\n", error.message[0] != '\0' ? error.message : tm_status_string(status));
        tm_domain_free(&domain);
        return 1;
    }

    status = tm_analyze_mesh(&mesh, &summary, &error);
    if (status != TM_STATUS_OK) {
        fprintf(stderr, "analyze failed: %s\n", error.message[0] != '\0' ? error.message : tm_status_string(status));
        tm_mesh_free(&mesh);
        tm_domain_free(&domain);
        return 1;
    }

    printf("triangle_count=%zu\n", summary.triangle_count);
    printf("point_count=%zu\n", summary.point_count);
    printf("input_point_count=%zu\n", summary.input_point_count);
    printf("steiner_point_count=%zu\n", summary.steiner_point_count);
    printf("area_min=%.17g\n", summary.area_min);
    printf("area_max=%.17g\n", summary.area_max);
    printf("min_angle_deg_min=%.17g\n", summary.min_angle_deg_min);
    printf("count_min_angle_lt_20=%zu\n", summary.count_min_angle_lt_20);

    tm_mesh_free(&mesh);
    tm_domain_free(&domain);
    return 0;
}
