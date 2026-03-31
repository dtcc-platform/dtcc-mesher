#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mesh.h"
#include "predicates.h"
#include "report.h"
#include "validate.h"
#include "io_pslg.h"
#include "cdt.h"

typedef struct {
    const char *name;
    size_t expected_triangle_count;
} SuccessCase;

static int verify_ccw(const TMMesh *mesh);
static int verify_delaunay(const TMMesh *mesh);
static int smoke_write_reports(
    const char *name,
    const TMMesh *mesh,
    const TMTriangleMetric *metrics,
    const TMSummary *summary
);

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
        fprintf(stderr, "missing output file %s\n", path);
        return 0;
    }
    if (st.st_size <= 0) {
        fprintf(stderr, "empty output file %s\n", path);
        return 0;
    }
    return 1;
}

static int build_case_path(char *buffer, size_t buffer_size, const char *name, const char *suffix)
{
    int needed = snprintf(buffer, buffer_size, "tests/%s/%s%s", "cases", name, suffix);
    return needed >= 0 && (size_t) needed < buffer_size;
}

static int build_tmp_path(char *buffer, size_t buffer_size, const char *name, const char *suffix)
{
    int needed = snprintf(buffer, buffer_size, "tests/tmp/%s%s", name, suffix);
    return needed >= 0 && (size_t) needed < buffer_size;
}

static int nearly_equal(double lhs, double rhs, double tolerance)
{
    return fabs(lhs - rhs) <= tolerance;
}

static int validate_mesh(const TMMesh *mesh, int check_delaunay)
{
    TMValidationReport report;
    TMStatus status = tm_validate_mesh(mesh, check_delaunay, &report);

    if (status != TM_OK) {
        fprintf(
            stderr,
            "mesh validation failed: adjacency=%zu orientation=%zu duplicates=%zu incident=%zu constrained=%zu delaunay=%zu\n",
            report.adjacency_errors,
            report.orientation_errors,
            report.duplicate_triangle_errors,
            report.incident_triangle_errors,
            report.constrained_edge_errors,
            report.local_delaunay_errors
        );
        return 0;
    }

    return 1;
}

static int validate_quality_mesh(const TMMesh *mesh, double min_angle_deg)
{
    TMValidationReport report;
    TMStatus status = tm_validate_quality_mesh(mesh, min_angle_deg, &report);

    if (status != TM_OK) {
        fprintf(
            stderr,
            "quality validation failed: adjacency=%zu orientation=%zu duplicates=%zu incident=%zu constrained=%zu delaunay=%zu encroached=%zu bad=%zu exempt=%zu\n",
            report.adjacency_errors,
            report.orientation_errors,
            report.duplicate_triangle_errors,
            report.incident_triangle_errors,
            report.constrained_edge_errors,
            report.local_delaunay_errors,
            report.encroached_segment_errors,
            report.bad_triangle_errors,
            report.exempt_triangle_count
        );
        return 0;
    }

    return 1;
}

static int load_case_mesh(const char *name, TMMesh *mesh)
{
    char case_path[256];
    TMPoint *points = NULL;
    size_t point_count = 0;
    TMStatus status;
    int ok = 0;

    memset(mesh, 0, sizeof(*mesh));

    if (!build_case_path(case_path, sizeof(case_path), name, ".pts")) {
        fprintf(stderr, "case path too long for %s\n", name);
        return 0;
    }

    status = tm_read_points_file(case_path, &points, &point_count);
    if (status != TM_OK) {
        fprintf(stderr, "failed to read %s: %s\n", name, tm_internal_status_string(status));
        goto cleanup;
    }

    status = tm_build_mesh(points, point_count, mesh);
    if (status != TM_OK) {
        fprintf(stderr, "failed to build mesh for %s: %s\n", name, tm_internal_status_string(status));
        goto cleanup;
    }

    ok = 1;

cleanup:
    tm_free_points(points);
    if (!ok) {
        tm_free_mesh(mesh);
    }
    return ok;
}

static int load_pslg_mesh_with_options(
    const char *name,
    int refine,
    int use_offcenters,
    int protect_acute_corners,
    TMAcuteProtectionMode acute_mode,
    double min_angle_deg,
    TMMesh *mesh
)
{
    char case_path[256];
    TMPSLG pslg;
    TMBuildOptions options;
    TMStatus status;
    int ok = 0;

    memset(mesh, 0, sizeof(*mesh));
    memset(&pslg, 0, sizeof(pslg));
    options.verbose = 0;
    options.refine = refine;
    options.use_offcenters = use_offcenters;
    options.protect_acute_corners = protect_acute_corners;
    options.acute_mode = acute_mode;
    options.min_angle_deg = min_angle_deg;
    options.protect_angle_deg = 0.0;
    options.max_refinement_steps = 0;
    options.max_protection_levels = 6;

    if (!build_case_path(case_path, sizeof(case_path), name, ".pslg")) {
        fprintf(stderr, "case path too long for %s\n", name);
        return 0;
    }

    status = tm_read_pslg_file(case_path, &pslg);
    if (status != TM_OK) {
        fprintf(stderr, "failed to read pslg %s: %s\n", name, tm_internal_status_string(status));
        goto cleanup;
    }

    status = tm_build_pslg_mesh(&pslg, &options, mesh);
    if (status != TM_OK) {
        fprintf(stderr, "failed to build pslg mesh for %s: %s\n", name, tm_internal_status_string(status));
        goto cleanup;
    }

    ok = 1;

cleanup:
    tm_free_pslg(&pslg);
    if (!ok) {
        tm_free_mesh(mesh);
    }
    return ok;
}

static int load_pslg_mesh(const char *name, TMMesh *mesh)
{
    return load_pslg_mesh_with_options(name, 0, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, mesh);
}

static int find_first_interior_edge(const TMMesh *mesh, int *out_triangle, int *out_edge)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge;

        for (edge = 0; edge < 3; ++edge) {
            if (mesh->triangles[tri_index].nbr[edge] > (int) tri_index) {
                *out_triangle = (int) tri_index;
                *out_edge = edge;
                return 1;
            }
        }
    }

    return 0;
}

static int point_in_triangle(const TMMesh *mesh, int triangle_index, const double point[2])
{
    int edge;

    for (edge = 0; edge < 3; ++edge) {
        int a;
        int b;
        double side;

        tm_triangle_edge_vertices(&mesh->triangles[triangle_index], edge, &a, &b);
        side = orient2d((REAL *) mesh->points[a].xy, (REAL *) mesh->points[b].xy, (REAL *) point);
        if (side < 0.0) {
            return 0;
        }
    }

    return 1;
}

static int edge_pair_exists(const TMMesh *mesh, int a, int b)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        if (tm_find_edge_in_triangle(&mesh->triangles[tri_index], a, b) >= 0) {
            return 1;
        }
    }

    return 0;
}

static int run_pslg_success_case(const char *name, size_t expected_triangle_count)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh(name, &mesh);

    if (!ok) {
        return 0;
    }

    ok = mesh.triangle_count == expected_triangle_count &&
         validate_mesh(&mesh, 1) &&
         verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "pslg success case %s failed basic checks\n", name);
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for pslg %s\n", name);
        ok = 0;
        goto cleanup;
    }

    if (!smoke_write_reports(name, &mesh, metrics, &summary)) {
        ok = 0;
    }

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int run_refined_pslg_quality_case(const char *name, double min_angle_deg, int require_steiner_points)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options(name, 1, 0, 1, TM_ACUTE_MODE_SHELL, min_angle_deg, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, min_angle_deg) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "refined PSLG case %s failed quality checks\n", name);
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for refined pslg %s\n", name);
        ok = 0;
        goto cleanup;
    }

    if (require_steiner_points && summary.steiner_point_count == 0) {
        fprintf(stderr, "refined PSLG case %s did not add Steiner points\n", name);
        ok = 0;
    }
    if (summary.count_min_angle_lt_20 != 0) {
        fprintf(stderr, "refined PSLG case %s still has min angles below 20 degrees\n", name);
        ok = 0;
    }
    if (!smoke_write_reports(name, &mesh, metrics, &summary)) {
        ok = 0;
    }

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int run_protected_pslg_quality_case(
    const char *name,
    double min_angle_deg,
    size_t min_protected_corners,
    size_t min_exempt_triangles
)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options(name, 1, 0, 1, TM_ACUTE_MODE_SHELL, min_angle_deg, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, min_angle_deg) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "protected PSLG case %s failed quality checks\n", name);
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for protected pslg %s\n", name);
        ok = 0;
        goto cleanup;
    }

    if (summary.protected_corner_count < min_protected_corners) {
        fprintf(stderr, "protected PSLG case %s did not register enough protected corners\n", name);
        ok = 0;
    }
    if (summary.exempt_triangle_count < min_exempt_triangles) {
        fprintf(stderr, "protected PSLG case %s did not create exempt shield triangles\n", name);
        ok = 0;
    }
    if (min_exempt_triangles != 0 && summary.count_min_angle_lt_20 == 0) {
        fprintf(stderr, "protected PSLG case %s should retain exempt low-angle triangles near acute corners\n", name);
        ok = 0;
    }
    if (!smoke_write_reports(name, &mesh, metrics, &summary)) {
        ok = 0;
    }

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int run_pslg_failure_case(const char *name, TMStatus expected_status)
{
    char case_path[256];
    TMPSLG pslg;
    TMStatus status;

    memset(&pslg, 0, sizeof(pslg));

    if (!build_case_path(case_path, sizeof(case_path), name, ".pslg")) {
        fprintf(stderr, "case path too long for %s\n", name);
        return 0;
    }

    status = tm_read_pslg_file(case_path, &pslg);
    tm_free_pslg(&pslg);

    if (status != expected_status) {
        fprintf(
            stderr,
            "unexpected PSLG status for %s: expected %s, got %s\n",
            name,
            tm_internal_status_string(expected_status),
            tm_internal_status_string(status)
        );
        return 0;
    }

    return 1;
}

static int verify_ccw(const TMMesh *mesh)
{
    size_t i;

    for (i = 0; i < mesh->triangle_count; ++i) {
        const TMTriangle *triangle = &mesh->triangles[i];
        double orientation = orient2d(
            (REAL *) mesh->points[triangle->v[0]].xy,
            (REAL *) mesh->points[triangle->v[1]].xy,
            (REAL *) mesh->points[triangle->v[2]].xy
        );

        if (orientation <= 0.0) {
            fprintf(stderr, "triangle %zu is not CCW\n", i);
            return 0;
        }
    }

    return 1;
}

static int verify_delaunay(const TMMesh *mesh)
{
    size_t tri_index;
    size_t point_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        const TMTriangle *triangle = &mesh->triangles[tri_index];

        for (point_index = 0; point_index < mesh->point_count; ++point_index) {
            double value;

            if ((int) point_index == triangle->v[0] ||
                (int) point_index == triangle->v[1] ||
                (int) point_index == triangle->v[2]) {
                continue;
            }

            value = incircle(
                (REAL *) mesh->points[triangle->v[0]].xy,
                (REAL *) mesh->points[triangle->v[1]].xy,
                (REAL *) mesh->points[triangle->v[2]].xy,
                (REAL *) mesh->points[point_index].xy
            );

            if (value > 0.0) {
                fprintf(stderr, "point %zu violates Delaunay condition for triangle %zu\n", point_index, tri_index);
                return 0;
            }
        }
    }

    return 1;
}

static int smoke_write_reports(
    const char *name,
    const TMMesh *mesh,
    const TMTriangleMetric *metrics,
    const TMSummary *summary
)
{
    char tri_path[256];
    char svg_path[256];
    char metrics_path[256];
    char summary_path[256];
    TMStatus status;

    if (!ensure_directory("tests/tmp")) {
        return 0;
    }

    if (!build_tmp_path(tri_path, sizeof(tri_path), name, ".tri") ||
        !build_tmp_path(svg_path, sizeof(svg_path), name, ".svg") ||
        !build_tmp_path(metrics_path, sizeof(metrics_path), name, ".metrics.csv") ||
        !build_tmp_path(summary_path, sizeof(summary_path), name, ".summary.txt")) {
        fprintf(stderr, "failed to build temp output paths for %s\n", name);
        return 0;
    }

    status = tm_write_tri_file(tri_path, mesh);
    if (status != TM_OK) {
        fprintf(stderr, "failed to write tri report for %s: %s\n", name, tm_internal_status_string(status));
        return 0;
    }
    status = tm_write_svg_file(svg_path, mesh);
    if (status != TM_OK) {
        fprintf(stderr, "failed to write svg report for %s: %s\n", name, tm_internal_status_string(status));
        return 0;
    }
    status = tm_write_metrics_csv(metrics_path, metrics, mesh->triangle_count);
    if (status != TM_OK) {
        fprintf(stderr, "failed to write metrics report for %s: %s\n", name, tm_internal_status_string(status));
        return 0;
    }
    status = tm_write_summary_file(summary_path, summary);
    if (status != TM_OK) {
        fprintf(stderr, "failed to write summary report for %s: %s\n", name, tm_internal_status_string(status));
        return 0;
    }

    return file_exists_and_nonempty(tri_path) &&
           file_exists_and_nonempty(svg_path) &&
           file_exists_and_nonempty(metrics_path) &&
           file_exists_and_nonempty(summary_path);
}

static int run_success_case(
    const SuccessCase *test_case,
    TMMesh *out_mesh,
    TMTriangleMetric **out_metrics,
    TMSummary *out_summary
)
{
    TMStatus status;
    int ok = 0;

    *out_metrics = NULL;

    if (!load_case_mesh(test_case->name, out_mesh)) {
        goto cleanup;
    }

    if (out_mesh->triangle_count != test_case->expected_triangle_count) {
        fprintf(
            stderr,
            "unexpected triangle count for %s: expected %zu, got %zu\n",
            test_case->name,
            test_case->expected_triangle_count,
            out_mesh->triangle_count
        );
        goto cleanup;
    }

    if (!validate_mesh(out_mesh, 1)) {
        goto cleanup;
    }

    if (!verify_ccw(out_mesh) || !verify_delaunay(out_mesh)) {
        goto cleanup;
    }

    status = tm_compute_metrics(out_mesh, out_metrics, out_summary);
    if (status != TM_OK) {
        fprintf(stderr, "failed to compute metrics for %s: %s\n", test_case->name, tm_internal_status_string(status));
        goto cleanup;
    }

    if (!smoke_write_reports(test_case->name, out_mesh, *out_metrics, out_summary)) {
        goto cleanup;
    }

    ok = 1;

cleanup:
    if (!ok) {
        tm_free_metrics(*out_metrics);
        *out_metrics = NULL;
        tm_free_mesh(out_mesh);
    }
    return ok;
}

static int run_failure_case(const char *name, TMStatus expected_status)
{
    char case_path[256];
    TMPoint *points = NULL;
    size_t point_count = 0;
    TMMesh mesh;
    TMStatus status;

    memset(&mesh, 0, sizeof(mesh));

    if (!build_case_path(case_path, sizeof(case_path), name, ".pts")) {
        fprintf(stderr, "case path too long for %s\n", name);
        return 0;
    }

    status = tm_read_points_file(case_path, &points, &point_count);
    if (status != TM_OK) {
        fprintf(stderr, "failed to read %s: %s\n", name, tm_internal_status_string(status));
        tm_free_points(points);
        return 0;
    }

    status = tm_build_mesh(points, point_count, &mesh);
    tm_free_points(points);
    tm_free_mesh(&mesh);

    if (status != expected_status) {
        fprintf(
            stderr,
            "unexpected status for %s: expected %s, got %s\n",
            name,
            tm_internal_status_string(expected_status),
            tm_internal_status_string(status)
        );
        return 0;
    }

    return 1;
}

static int test_triangle3(void)
{
    SuccessCase test_case = {"triangle3", 1};
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok;

    ok = run_success_case(&test_case, &mesh, &metrics, &summary);
    if (!ok) {
        return 0;
    }

    ok = nearly_equal(metrics[0].min_angle_deg, 45.0, 1e-6) &&
         nearly_equal(metrics[0].max_angle_deg, 90.0, 1e-6);
    if (!ok) {
        fprintf(stderr, "triangle3 diagnostics are off\n");
    }

    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_square4(void)
{
    SuccessCase test_case = {"square4", 2};
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = run_success_case(&test_case, &mesh, &metrics, &summary);

    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_square_center5(void)
{
    SuccessCase test_case = {"square_center5", 4};
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    size_t i;
    int ok = run_success_case(&test_case, &mesh, &metrics, &summary);

    if (!ok) {
        return 0;
    }

    for (i = 0; i < mesh.triangle_count; ++i) {
        const TMTriangle *triangle = &mesh.triangles[i];

        if (triangle->v[0] != 4 && triangle->v[1] != 4 && triangle->v[2] != 4) {
            fprintf(stderr, "square_center5 triangle %zu does not use center point\n", i);
            ok = 0;
            break;
        }
    }

    if (ok && !nearly_equal(summary.min_angle_deg_min, 45.0, 1e-6)) {
        fprintf(stderr, "square_center5 global min angle is off\n");
        ok = 0;
    }

    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_thin_rect4(void)
{
    SuccessCase test_case = {"thin_rect4", 2};
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = run_success_case(&test_case, &mesh, &metrics, &summary);

    if (ok && !(summary.min_angle_deg_min < 10.0)) {
        fprintf(stderr, "thin_rect4 min angle should be below 10 degrees\n");
        ok = 0;
    }

    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_warped_disc64(void)
{
    SuccessCase test_case = {"warped_disc64", 112};
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = run_success_case(&test_case, &mesh, &metrics, &summary);

    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_point_location_walk(void)
{
    TMMesh mesh;
    TMLocation location;
    double query[2] = {0.2, 0.1};
    TMStatus status;
    int ok = load_case_mesh("square_center5", &mesh);

    if (!ok) {
        return 0;
    }

    status = tm_locate_point(&mesh, query, (int) mesh.triangle_count - 1, &location);
    ok = status == TM_OK &&
         location.triangle >= 0 &&
         (size_t) location.triangle < mesh.triangle_count &&
         point_in_triangle(&mesh, location.triangle, query);
    if (!ok) {
        fprintf(stderr, "point location walk failed\n");
    }

    tm_free_mesh(&mesh);
    return ok;
}

static int test_edge_flip_primitive(void)
{
    TMMesh mesh;
    int triangle_index;
    int edge;
    int a;
    int b;
    TMStatus status;
    int ok = load_case_mesh("square4", &mesh);

    if (!ok) {
        return 0;
    }

    ok = find_first_interior_edge(&mesh, &triangle_index, &edge);
    if (!ok) {
        fprintf(stderr, "failed to find an interior edge to flip\n");
        tm_free_mesh(&mesh);
        return 0;
    }

    tm_triangle_edge_vertices(&mesh.triangles[triangle_index], edge, &a, &b);
    status = tm_flip_edge(&mesh, triangle_index, edge);
    ok = status == TM_OK &&
         mesh.triangle_count == 2 &&
         !edge_pair_exists(&mesh, a, b) &&
         validate_mesh(&mesh, 1) &&
         verify_ccw(&mesh) &&
         verify_delaunay(&mesh);
    if (!ok) {
        fprintf(stderr, "edge flip primitive failed\n");
    }

    tm_free_mesh(&mesh);
    return ok;
}

static int test_insert_point_in_triangle_primitive(void)
{
    TMMesh mesh;
    double query[2] = {0.25, 0.25};
    int point_index = -1;
    TMStatus status;
    int ok = load_case_mesh("triangle3", &mesh);

    if (!ok) {
        return 0;
    }

    status = tm_insert_point_in_triangle(&mesh, 0, query, TM_VERTEX_TRIANGLE_SPLIT, &point_index);
    ok = status == TM_OK &&
         mesh.point_count == 4 &&
         mesh.triangle_count == 3 &&
         point_index >= 0 &&
         mesh.points[point_index].kind == TM_VERTEX_TRIANGLE_SPLIT &&
         validate_mesh(&mesh, 1) &&
         verify_ccw(&mesh) &&
         verify_delaunay(&mesh);
    if (!ok) {
        fprintf(stderr, "triangle insertion primitive failed\n");
    }

    tm_free_mesh(&mesh);
    return ok;
}

static int test_insert_point_on_edge_primitive(void)
{
    TMMesh mesh;
    int triangle_index;
    int edge;
    int a;
    int b;
    double query[2];
    int point_index = -1;
    TMStatus status;
    int ok = load_case_mesh("square4", &mesh);

    if (!ok) {
        return 0;
    }

    ok = find_first_interior_edge(&mesh, &triangle_index, &edge);
    if (!ok) {
        fprintf(stderr, "failed to find an interior edge to split\n");
        tm_free_mesh(&mesh);
        return 0;
    }

    tm_triangle_edge_vertices(&mesh.triangles[triangle_index], edge, &a, &b);
    query[0] = 0.5 * (mesh.points[a].xy[0] + mesh.points[b].xy[0]);
    query[1] = 0.5 * (mesh.points[a].xy[1] + mesh.points[b].xy[1]);

    status = tm_insert_point_on_edge(&mesh, triangle_index, edge, query, TM_VERTEX_SEGMENT_SPLIT, &point_index);
    ok = status == TM_OK &&
         mesh.point_count == 5 &&
         mesh.triangle_count == 4 &&
         point_index >= 0 &&
         mesh.points[point_index].kind == TM_VERTEX_SEGMENT_SPLIT &&
         validate_mesh(&mesh, 1) &&
         verify_ccw(&mesh) &&
         verify_delaunay(&mesh);
    if (!ok) {
        fprintf(stderr, "edge insertion primitive failed\n");
    }

    tm_free_mesh(&mesh);
    return ok;
}

static int test_square_pslg(void)
{
    return run_pslg_success_case("square_domain", 2);
}

static int test_l_shape_pslg(void)
{
    return run_pslg_success_case("l_shape_domain", 4);
}

static int test_square_hole_pslg(void)
{
    return run_pslg_success_case("square_hole_domain", 8);
}

static int test_multi_hole_pslg(void)
{
    return run_pslg_success_case("multi_hole_domain", 14);
}

static int test_warped_disc64_domain_pslg(void)
{
    return run_pslg_success_case("warped_disc64_domain", 108);
}

static int test_square_hole_pslg_refined(void)
{
    return run_refined_pslg_quality_case("square_hole_domain", 20.0, 1);
}

static int test_long_channel_pslg_refined(void)
{
    return run_refined_pslg_quality_case("long_channel_domain", 20.0, 1);
}

static int test_top_cluster_pslg_refined(void)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options("top_cluster_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, 20.0) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "refined PSLG case top_cluster_domain failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for refined pslg top_cluster_domain\n");
        ok = 0;
        goto cleanup;
    }

    if (summary.triangle_split_point_count == 0) {
        fprintf(stderr, "refined PSLG case top_cluster_domain did not exercise triangle split insertions\n");
        ok = 0;
    }
    if (!smoke_write_reports("top_cluster_domain", &mesh, metrics, &summary)) {
        ok = 0;
    }

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_top_cluster_pslg_offcenters(void)
{
    TMMesh circum_mesh;
    TMMesh offcenter_mesh;
    TMTriangleMetric *circum_metrics = NULL;
    TMTriangleMetric *offcenter_metrics = NULL;
    TMSummary circum_summary;
    TMSummary offcenter_summary;
    int ok = 0;

    memset(&circum_mesh, 0, sizeof(circum_mesh));
    memset(&offcenter_mesh, 0, sizeof(offcenter_mesh));

    ok = load_pslg_mesh_with_options("top_cluster_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &circum_mesh) &&
         load_pslg_mesh_with_options("top_cluster_domain", 1, 1, 1, TM_ACUTE_MODE_SHELL, 20.0, &offcenter_mesh);
    if (!ok) {
        tm_free_mesh(&circum_mesh);
        tm_free_mesh(&offcenter_mesh);
        return 0;
    }

    ok = validate_quality_mesh(&circum_mesh, 20.0) &&
         validate_quality_mesh(&offcenter_mesh, 20.0) &&
         verify_ccw(&circum_mesh) &&
         verify_ccw(&offcenter_mesh);
    if (!ok) {
        fprintf(stderr, "top_cluster_domain off-center comparison failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&circum_mesh, &circum_metrics, &circum_summary) != TM_OK ||
        tm_compute_metrics(&offcenter_mesh, &offcenter_metrics, &offcenter_summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for off-center comparison\n");
        ok = 0;
        goto cleanup;
    }

    ok = offcenter_summary.steiner_point_count <= circum_summary.steiner_point_count &&
         offcenter_summary.triangle_count <= circum_summary.triangle_count &&
         offcenter_summary.count_min_angle_lt_20 == 0;
    if (!ok) {
        fprintf(stderr, "off-center mode did not improve top_cluster_domain as expected\n");
        goto cleanup;
    }

    ok = smoke_write_reports("top_cluster_domain_offcenter", &offcenter_mesh, offcenter_metrics, &offcenter_summary);

cleanup:
    tm_free_metrics(circum_metrics);
    tm_free_metrics(offcenter_metrics);
    tm_free_mesh(&circum_mesh);
    tm_free_mesh(&offcenter_mesh);
    return ok;
}

static int test_warped_disc64_domain_pslg_refined(void)
{
    return run_refined_pslg_quality_case("warped_disc64_domain", 20.0, 1);
}

static int test_city_footprints_domain_refined(void)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options("city_footprints_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, 20.0) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "city_footprints_domain failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for city_footprints_domain\n");
        ok = 0;
        goto cleanup;
    }

    ok = summary.input_point_count >= 180 &&
         summary.steiner_point_count >= 200 &&
         summary.segment_split_point_count >= 80 &&
         summary.triangle_split_point_count >= 120 &&
         summary.triangle_count >= 550 &&
         summary.area_min <= 1.0 &&
         summary.area_max >= 200.0 &&
         summary.area_max >= summary.area_min * 200.0;
    if (!ok) {
        fprintf(stderr, "city_footprints_domain did not show the expected grading stress pattern\n");
        goto cleanup;
    }

    ok = smoke_write_reports("city_footprints_domain", &mesh, metrics, &summary);

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_city_downtown_domain_refined(void)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options("city_downtown_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, 20.0) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "city_downtown_domain failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for city_downtown_domain\n");
        ok = 0;
        goto cleanup;
    }

    ok = summary.input_point_count >= 180 &&
         summary.steiner_point_count >= 220 &&
         summary.segment_split_point_count >= 60 &&
         summary.triangle_split_point_count >= 150 &&
         summary.triangle_count >= 620 &&
         summary.area_min <= 1.0 &&
         summary.area_max >= 200.0 &&
         summary.area_max >= summary.area_min * 200.0;
    if (!ok) {
        fprintf(stderr, "city_downtown_domain did not show the expected realistic grading stress pattern\n");
        goto cleanup;
    }

    ok = smoke_write_reports("city_downtown_domain", &mesh, metrics, &summary);

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_city_tight_downtown_domain_refined(void)
{
    TMMesh mesh;
    TMTriangleMetric *metrics = NULL;
    TMSummary summary;
    int ok = load_pslg_mesh_with_options("city_tight_downtown_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &mesh);

    if (!ok) {
        return 0;
    }

    ok = validate_quality_mesh(&mesh, 20.0) && verify_ccw(&mesh);
    if (!ok) {
        fprintf(stderr, "city_tight_downtown_domain failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&mesh, &metrics, &summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for city_tight_downtown_domain\n");
        ok = 0;
        goto cleanup;
    }

    ok = summary.input_point_count >= 290 &&
         summary.steiner_point_count >= 450 &&
         summary.segment_split_point_count >= 300 &&
         summary.triangle_split_point_count >= 120 &&
         summary.triangle_count >= 950 &&
         summary.area_min <= 0.3 &&
         summary.area_max >= 1500.0 &&
         summary.area_max >= summary.area_min * 1000.0;
    if (!ok) {
        fprintf(stderr, "city_tight_downtown_domain did not show the expected tiny-gap grading stress pattern\n");
        goto cleanup;
    }

    ok = smoke_write_reports("city_tight_downtown_domain", &mesh, metrics, &summary);

cleanup:
    tm_free_metrics(metrics);
    tm_free_mesh(&mesh);
    return ok;
}

static int test_tiny_slit_shell_vs_simple(void)
{
    TMMesh simple_mesh;
    TMMesh shell_mesh;
    TMTriangleMetric *simple_metrics = NULL;
    TMTriangleMetric *shell_metrics = NULL;
    TMSummary simple_summary;
    TMSummary shell_summary;
    int ok = 0;

    memset(&simple_mesh, 0, sizeof(simple_mesh));
    memset(&shell_mesh, 0, sizeof(shell_mesh));

    ok = load_pslg_mesh_with_options("tiny_slit_domain", 1, 0, 1, TM_ACUTE_MODE_SIMPLE, 20.0, &simple_mesh) &&
         load_pslg_mesh_with_options("tiny_slit_domain", 1, 0, 1, TM_ACUTE_MODE_SHELL, 20.0, &shell_mesh);
    if (!ok) {
        tm_free_mesh(&simple_mesh);
        tm_free_mesh(&shell_mesh);
        return 0;
    }

    ok = validate_quality_mesh(&simple_mesh, 20.0) &&
         validate_quality_mesh(&shell_mesh, 20.0) &&
         verify_ccw(&simple_mesh) &&
         verify_ccw(&shell_mesh);
    if (!ok) {
        fprintf(stderr, "tiny_slit_domain simple vs shell comparison failed quality checks\n");
        goto cleanup;
    }

    if (tm_compute_metrics(&simple_mesh, &simple_metrics, &simple_summary) != TM_OK ||
        tm_compute_metrics(&shell_mesh, &shell_metrics, &shell_summary) != TM_OK) {
        fprintf(stderr, "failed to compute metrics for tiny_slit_domain comparison\n");
        ok = 0;
        goto cleanup;
    }

    ok = shell_summary.steiner_point_count < simple_summary.steiner_point_count &&
         shell_summary.triangle_count < simple_summary.triangle_count &&
         shell_summary.count_min_angle_lt_20 == 0;
    if (!ok) {
        fprintf(stderr, "shell acute mode did not improve tiny_slit_domain as expected\n");
        goto cleanup;
    }

    ok = smoke_write_reports("tiny_slit_domain_simple", &simple_mesh, simple_metrics, &simple_summary) &&
         smoke_write_reports("tiny_slit_domain_shell", &shell_mesh, shell_metrics, &shell_summary);

cleanup:
    tm_free_metrics(simple_metrics);
    tm_free_metrics(shell_metrics);
    tm_free_mesh(&simple_mesh);
    tm_free_mesh(&shell_mesh);
    return ok;
}

static int test_acute_wedge_protected(void)
{
    return run_protected_pslg_quality_case("acute_wedge_5deg", 20.0, 1, 1);
}

static int test_narrow_notch_protected(void)
{
    return run_protected_pslg_quality_case("narrow_notch_domain", 20.0, 1, 0);
}

static int test_double_notch_protected(void)
{
    return run_protected_pslg_quality_case("double_notch_domain", 20.0, 2, 0);
}

static int test_acute_wedge_needs_protection(void)
{
    char case_path[256];
    TMPSLG pslg;
    TMMesh mesh;
    TMBuildOptions options;
    TMStatus status;
    int ok = 0;

    memset(&pslg, 0, sizeof(pslg));
    memset(&mesh, 0, sizeof(mesh));
    options.verbose = 0;
    options.refine = 1;
    options.use_offcenters = 0;
    options.protect_acute_corners = 0;
    options.acute_mode = TM_ACUTE_MODE_SHELL;
    options.min_angle_deg = 20.0;
    options.protect_angle_deg = 0.0;
    options.max_refinement_steps = 0;
    options.max_protection_levels = 6;

    if (!build_case_path(case_path, sizeof(case_path), "acute_wedge_5deg", ".pslg")) {
        fprintf(stderr, "case path too long for acute_wedge_5deg\n");
        return 0;
    }

    status = tm_read_pslg_file(case_path, &pslg);
    if (status != TM_OK) {
        fprintf(stderr, "failed to read acute_wedge_5deg for negative protection test: %s\n", tm_internal_status_string(status));
        goto cleanup;
    }

    status = tm_build_pslg_mesh(&pslg, &options, &mesh);
    ok = status != TM_OK;
    if (!ok) {
        fprintf(stderr, "acute_wedge_5deg unexpectedly succeeded without acute protection\n");
    }

cleanup:
    tm_free_pslg(&pslg);
    tm_free_mesh(&mesh);
    return ok;
}

int main(void)
{
    int failures = 0;

    tm_initialize();

    if (!test_triangle3()) {
        failures += 1;
    }
    if (!test_square4()) {
        failures += 1;
    }
    if (!test_square_center5()) {
        failures += 1;
    }
    if (!test_thin_rect4()) {
        failures += 1;
    }
    if (!test_warped_disc64()) {
        failures += 1;
    }
    if (!test_point_location_walk()) {
        failures += 1;
    }
    if (!test_edge_flip_primitive()) {
        failures += 1;
    }
    if (!test_insert_point_in_triangle_primitive()) {
        failures += 1;
    }
    if (!test_insert_point_on_edge_primitive()) {
        failures += 1;
    }
    if (!test_square_pslg()) {
        failures += 1;
    }
    if (!test_l_shape_pslg()) {
        failures += 1;
    }
    if (!test_square_hole_pslg()) {
        failures += 1;
    }
    if (!test_multi_hole_pslg()) {
        failures += 1;
    }
    if (!test_warped_disc64_domain_pslg()) {
        failures += 1;
    }
    if (!test_square_hole_pslg_refined()) {
        failures += 1;
    }
    if (!test_long_channel_pslg_refined()) {
        failures += 1;
    }
    if (!test_top_cluster_pslg_refined()) {
        failures += 1;
    }
    if (!test_top_cluster_pslg_offcenters()) {
        failures += 1;
    }
    if (!test_warped_disc64_domain_pslg_refined()) {
        failures += 1;
    }
    if (!test_city_footprints_domain_refined()) {
        failures += 1;
    }
    if (!test_city_downtown_domain_refined()) {
        failures += 1;
    }
    if (!test_city_tight_downtown_domain_refined()) {
        failures += 1;
    }
    if (!test_tiny_slit_shell_vs_simple()) {
        failures += 1;
    }
    if (!test_acute_wedge_protected()) {
        failures += 1;
    }
    if (!test_narrow_notch_protected()) {
        failures += 1;
    }
    if (!test_double_notch_protected()) {
        failures += 1;
    }
    if (!test_acute_wedge_needs_protection()) {
        failures += 1;
    }
    if (!run_failure_case("duplicate5", TM_ERR_DUPLICATE)) {
        failures += 1;
    }
    if (!run_failure_case("collinear4", TM_ERR_COLLINEAR)) {
        failures += 1;
    }
    if (!run_pslg_failure_case("intersecting_segments", TM_ERR_INVALID_PSLG)) {
        failures += 1;
    }

    if (failures != 0) {
        fprintf(stderr, "%d test(s) failed\n", failures);
        return EXIT_FAILURE;
    }

    printf("all tests passed\n");
    return EXIT_SUCCESS;
}
