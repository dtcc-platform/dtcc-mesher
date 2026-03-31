#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"

static int ensure_directory(const char *path)
{
    if (mkdir(path, 0777) == 0 || errno == EEXIST) {
        return 1;
    }

    fprintf(stderr, "failed to create directory %s\n", path);
    return 0;
}

static int file_exists_and_nonempty(const char *path)
{
    struct stat st;

    if (stat(path, &st) != 0) {
        fprintf(stderr, "missing file %s\n", path);
        return 0;
    }

    if (st.st_size <= 0) {
        fprintf(stderr, "empty file %s\n", path);
        return 0;
    }

    return 1;
}

static int nearly_equal(double lhs, double rhs, double tolerance)
{
    return fabs(lhs - rhs) <= tolerance;
}

static int test_generate_point_mesh(void)
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
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_error error;
    dtcc_mesher_status status;

    memset(&mesh, 0, sizeof(mesh));
    memset(&error, 0, sizeof(error));
    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = NULL;
    domain.num_segments = 0;
    domain.holes = NULL;
    domain.num_holes = 0;

    dtcc_mesher_options_init(&options);
    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "point generate failed: %s\n", error.message);
        return 0;
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "point analyze failed: %s\n", error.message);
        dtcc_mesher_mesh_free(&mesh);
        return 0;
    }

    if (mesh.num_triangles != 4 || summary.count_min_angle_lt_20 != 0 ||
        !nearly_equal(summary.min_angle_deg_min, 45.0, 1e-6)) {
        fprintf(stderr, "unexpected point mesh summary\n");
        dtcc_mesher_mesh_free(&mesh);
        return 0;
    }

    dtcc_mesher_mesh_free(&mesh);
    return 1;
}

static int test_read_generate_and_write_pslg(void)
{
    char svg_path[512];
    char summary_path[512];
    char tri_path[512];
    char csv_path[512];
    dtcc_mesher_domain domain;
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_quality_summary summary;
    dtcc_mesher_error error;
    dtcc_mesher_status status;

    memset(&domain, 0, sizeof(domain));
    memset(&mesh, 0, sizeof(mesh));
    memset(&error, 0, sizeof(error));

    if (!ensure_directory("tests/tmp")) {
        return 0;
    }

    status = dtcc_mesher_read_domain_file(DTCC_MESHER_SOURCE_DIR "/tests/cases/square_hole_domain.pslg", &domain, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "read pslg failed: %s\n", error.message);
        return 0;
    }

    dtcc_mesher_options_init(&options);
    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "pslg generate failed: %s\n", error.message);
        dtcc_mesher_domain_free(&domain);
        return 0;
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "pslg analyze failed: %s\n", error.message);
        dtcc_mesher_mesh_free(&mesh);
        dtcc_mesher_domain_free(&domain);
        return 0;
    }

    snprintf(svg_path, sizeof(svg_path), "tests/tmp/api_square_hole.svg");
    snprintf(summary_path, sizeof(summary_path), "tests/tmp/api_square_hole.summary.txt");
    snprintf(tri_path, sizeof(tri_path), "tests/tmp/api_square_hole.tri");
    snprintf(csv_path, sizeof(csv_path), "tests/tmp/api_square_hole.metrics.csv");

    if (dtcc_mesher_write_svg(&mesh, svg_path, &error) != DTCC_MESHER_STATUS_OK ||
        dtcc_mesher_write_quality_summary(&mesh, summary_path, &error) != DTCC_MESHER_STATUS_OK ||
        dtcc_mesher_write_triangles(&mesh, tri_path, &error) != DTCC_MESHER_STATUS_OK ||
        dtcc_mesher_write_quality_csv(&mesh, csv_path, &error) != DTCC_MESHER_STATUS_OK) {
        fprintf(stderr, "public writer failed: %s\n", error.message);
        dtcc_mesher_mesh_free(&mesh);
        dtcc_mesher_domain_free(&domain);
        return 0;
    }

    if (!file_exists_and_nonempty(svg_path) || !file_exists_and_nonempty(summary_path) ||
        !file_exists_and_nonempty(tri_path) || !file_exists_and_nonempty(csv_path)) {
        dtcc_mesher_mesh_free(&mesh);
        dtcc_mesher_domain_free(&domain);
        return 0;
    }

    if (summary.triangle_count == 0 || summary.point_count < summary.input_point_count) {
        fprintf(stderr, "unexpected PSLG summary values\n");
        dtcc_mesher_mesh_free(&mesh);
        dtcc_mesher_domain_free(&domain);
        return 0;
    }

    dtcc_mesher_mesh_free(&mesh);
    dtcc_mesher_domain_free(&domain);
    return 1;
}

static int test_invalid_geometry_rejected(void)
{
    dtcc_mesher_point points[] = {
        {0.0, 0.0},
        {1.0, 0.0}
    };
    dtcc_mesher_domain domain;
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh;
    dtcc_mesher_error error;
    dtcc_mesher_status status;

    memset(&mesh, 0, sizeof(mesh));
    memset(&error, 0, sizeof(error));
    domain.points = points;
    domain.num_points = sizeof(points) / sizeof(points[0]);
    domain.segments = NULL;
    domain.num_segments = 0;
    domain.holes = NULL;
    domain.num_holes = 0;

    dtcc_mesher_options_init(&options);
    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    dtcc_mesher_mesh_free(&mesh);

    if (status != DTCC_MESHER_STATUS_GEOMETRY) {
        fprintf(stderr, "expected geometry rejection, got %s (%s)\n", dtcc_mesher_status_string(status), error.message);
        return 0;
    }

    return 1;
}

int main(void)
{
    int ok = 1;

    ok = test_generate_point_mesh() && ok;
    ok = test_read_generate_and_write_pslg() && ok;
    ok = test_invalid_geometry_rejected() && ok;

    if (!ok) {
        return 1;
    }

    printf("public C API tests passed\n");
    return 0;
}
