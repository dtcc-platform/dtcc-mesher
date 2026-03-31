#include <stdio.h>
#include <string.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"

int main(int argc, char **argv)
{
    dtcc_mesher_domain domain;
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_error error;
    dtcc_mesher_status status;

    if (argc != 2) {
        fprintf(stderr, "usage: api_dump input.(pts|pslg)\n");
        return 2;
    }

    memset(&domain, 0, sizeof(domain));
    memset(&mesh, 0, sizeof(mesh));
    memset(&summary, 0, sizeof(summary));
    memset(&error, 0, sizeof(error));

    status = dtcc_mesher_read_domain_file(argv[1], &domain, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "read failed: %s\n", error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
        return 1;
    }

    dtcc_mesher_options_init(&options);
    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "generate failed: %s\n", error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
        dtcc_mesher_domain_free(&domain);
        return 1;
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "analyze failed: %s\n", error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
        dtcc_mesher_mesh_free(&mesh);
        dtcc_mesher_domain_free(&domain);
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

    dtcc_mesher_mesh_free(&mesh);
    dtcc_mesher_domain_free(&domain);
    return 0;
}
