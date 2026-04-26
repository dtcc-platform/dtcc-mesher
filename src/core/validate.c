#include "validate.h"

#include <math.h>
#include <stdlib.h>

#include "predicates.h"

static const double tm_pi = 3.14159265358979323846;

static void tm_sort_three_ints(int values[3])
{
    int tmp;

    if (values[0] > values[1]) {
        tmp = values[0];
        values[0] = values[1];
        values[1] = tmp;
    }
    if (values[1] > values[2]) {
        tmp = values[1];
        values[1] = values[2];
        values[2] = tmp;
    }
    if (values[0] > values[1]) {
        tmp = values[0];
        values[0] = values[1];
        values[1] = tmp;
    }
}

TMStatus tm_validate_mesh(const TMMesh *mesh, int check_delaunay, TMValidationReport *out_report)
{
    size_t tri_index;
    unsigned char *used_points = NULL;
    TMStatus status = TM_OK;

    if (mesh == NULL || mesh->points == NULL || mesh->triangles == NULL || out_report == NULL) {
        return TM_ERR_INTERNAL;
    }

    used_points = (unsigned char *) calloc(mesh->point_count, sizeof(*used_points));
    if (used_points == NULL && mesh->point_count != 0) {
        return TM_ERR_ALLOC;
    }

    out_report->triangle_count = mesh->triangle_count;
    out_report->boundary_edge_count = 0;
    out_report->adjacency_errors = 0;
    out_report->orientation_errors = 0;
    out_report->duplicate_triangle_errors = 0;
    out_report->incident_triangle_errors = 0;
    out_report->constrained_edge_errors = 0;
    out_report->local_delaunay_errors = 0;
    out_report->encroached_segment_errors = 0;
    out_report->bad_triangle_errors = 0;
    out_report->exempt_triangle_count = 0;
    out_report->first_encroached_segment = (size_t) -1;
    out_report->first_encroaching_point = -1;
    out_report->first_bad_triangle = (size_t) -1;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        const TMTriangle *triangle = &mesh->triangles[tri_index];
        double orientation = orient2d(
            (REAL *) mesh->points[triangle->v[0]].xy,
            (REAL *) mesh->points[triangle->v[1]].xy,
            (REAL *) mesh->points[triangle->v[2]].xy
        );
        int edge;

        used_points[triangle->v[0]] = 1;
        used_points[triangle->v[1]] = 1;
        used_points[triangle->v[2]] = 1;

        if (orientation <= 0.0) {
            out_report->orientation_errors += 1;
        }

        for (edge = 0; edge < 3; ++edge) {
            int neighbor_index = triangle->nbr[edge];

            if (neighbor_index < 0) {
                out_report->boundary_edge_count += 1;
                continue;
            }

            if ((size_t) neighbor_index >= mesh->triangle_count) {
                out_report->adjacency_errors += 1;
                continue;
            }

            {
                int a;
                int b;
                int neighbor_edge;
                const TMTriangle *neighbor = &mesh->triangles[neighbor_index];

                tm_triangle_edge_vertices(triangle, edge, &a, &b);
                neighbor_edge = tm_find_edge_in_triangle(neighbor, a, b);
                if (neighbor_edge < 0 || neighbor->nbr[neighbor_edge] != (int) tri_index) {
                    out_report->adjacency_errors += 1;
                    continue;
                }

                if (triangle->constrained[edge] != neighbor->constrained[neighbor_edge]) {
                    out_report->constrained_edge_errors += 1;
                }

                if (check_delaunay && tri_index < (size_t) neighbor_index && !tm_edge_is_locally_delaunay(mesh, (int) tri_index, edge)) {
                    out_report->local_delaunay_errors += 1;
                }
            }
        }
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        size_t other_index;
        int tri_key[3] = {
            mesh->triangles[tri_index].v[0],
            mesh->triangles[tri_index].v[1],
            mesh->triangles[tri_index].v[2]
        };

        tm_sort_three_ints(tri_key);

        for (other_index = tri_index + 1; other_index < mesh->triangle_count; ++other_index) {
            int other_key[3] = {
                mesh->triangles[other_index].v[0],
                mesh->triangles[other_index].v[1],
                mesh->triangles[other_index].v[2]
            };

            tm_sort_three_ints(other_key);
            if (tri_key[0] == other_key[0] && tri_key[1] == other_key[1] && tri_key[2] == other_key[2]) {
                out_report->duplicate_triangle_errors += 1;
            }
        }
    }

    for (tri_index = 0; tri_index < mesh->point_count; ++tri_index) {
        int incident_triangle = mesh->points[tri_index].incident_triangle;

        if (!used_points[tri_index]) {
            if (incident_triangle != -1) {
                out_report->incident_triangle_errors += 1;
            }
            continue;
        }

        if (incident_triangle < 0 || (size_t) incident_triangle >= mesh->triangle_count) {
            out_report->incident_triangle_errors += 1;
            continue;
        }

        if (!tm_triangle_contains_vertex(&mesh->triangles[incident_triangle], (int) tri_index)) {
            out_report->incident_triangle_errors += 1;
        }
    }

    for (tri_index = 0; tri_index < mesh->segment_count; ++tri_index) {
        size_t find_index;
        int found = 0;

        if (!mesh->segments[tri_index].live) {
            continue;
        }

        for (find_index = 0; find_index < mesh->triangle_count; ++find_index) {
            int edge = tm_find_edge_in_triangle(
                &mesh->triangles[find_index],
                mesh->segments[tri_index].v[0],
                mesh->segments[tri_index].v[1]
            );

            if (edge >= 0) {
                if (!mesh->triangles[find_index].constrained[edge]) {
                    out_report->constrained_edge_errors += 1;
                }
                found = 1;
            }
        }

        if (!found) {
            out_report->constrained_edge_errors += 1;
        }
    }

    if (out_report->adjacency_errors != 0 ||
        out_report->orientation_errors != 0 ||
        out_report->duplicate_triangle_errors != 0 ||
        out_report->incident_triangle_errors != 0 ||
        out_report->constrained_edge_errors != 0 ||
        out_report->local_delaunay_errors != 0) {
        status = TM_ERR_INVALID_MESH;
    }

    free(used_points);
    return status;
}

static double tm_distance_squared(const double a[2], const double b[2])
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];

    return dx * dx + dy * dy;
}

static int tm_point_encroaches_segment(
    const TMMesh *mesh,
    const TMSegment *segment,
    int point_index
)
{
    const double *a = mesh->points[segment->v[0]].xy;
    const double *b = mesh->points[segment->v[1]].xy;
    const double *p;
    double length_sq;
    double tolerance;
    double dot;

    if (point_index < 0 || (size_t) point_index >= mesh->point_count) {
        return 0;
    }
    if (point_index == segment->v[0] || point_index == segment->v[1]) {
        return 0;
    }
    if (mesh->points[point_index].incident_triangle < 0) {
        return 0;
    }

    p = mesh->points[point_index].xy;
    length_sq = tm_distance_squared(a, b);
    tolerance = 1e-12 * ((length_sq > 1.0) ? length_sq : 1.0);
    dot = (p[0] - a[0]) * (p[0] - b[0]) + (p[1] - a[1]) * (p[1] - b[1]);
    return dot < -tolerance;
}

static int tm_segment_triangle_edge(
    const TMMesh *mesh,
    const TMSegment *segment,
    int *out_triangle_index,
    int *out_edge
)
{
    size_t tri_index;

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        int edge = tm_find_edge_in_triangle(
            &mesh->triangles[tri_index],
            segment->v[0],
            segment->v[1]
        );

        if (edge >= 0) {
            if (out_triangle_index != NULL) {
                *out_triangle_index = (int) tri_index;
            }
            if (out_edge != NULL) {
                *out_edge = edge;
            }
            return 1;
        }
    }

    return 0;
}

static int tm_segment_is_encroached(
    const TMMesh *mesh,
    size_t segment_index,
    int *out_point_index
)
{
    const TMSegment *segment = &mesh->segments[segment_index];
    int triangle_index;
    int edge;

    if (out_point_index != NULL) {
        *out_point_index = -1;
    }
    if (!tm_segment_triangle_edge(mesh, segment, &triangle_index, &edge)) {
        return 0;
    }

    /* Only vertices on triangulated side(s) are visible to the constrained segment. */
    if (tm_point_encroaches_segment(mesh, segment, mesh->triangles[triangle_index].v[edge])) {
        if (out_point_index != NULL) {
            *out_point_index = mesh->triangles[triangle_index].v[edge];
        }
        return 1;
    }

    {
        int neighbor_index = mesh->triangles[triangle_index].nbr[edge];

        if (neighbor_index >= 0 && (size_t) neighbor_index < mesh->triangle_count) {
            int neighbor_edge = tm_find_edge_in_triangle(
                &mesh->triangles[neighbor_index],
                segment->v[0],
                segment->v[1]
            );

            if (neighbor_edge >= 0 &&
                tm_point_encroaches_segment(mesh, segment, mesh->triangles[neighbor_index].v[neighbor_edge])) {
                if (out_point_index != NULL) {
                    *out_point_index = mesh->triangles[neighbor_index].v[neighbor_edge];
                }
                return 1;
            }
        }
    }

    return 0;
}

static int tm_triangle_is_bad(const TMMesh *mesh, size_t tri_index, double beta)
{
    const TMTriangle *triangle = &mesh->triangles[tri_index];
    double len01 = sqrt(tm_distance_squared(mesh->points[triangle->v[0]].xy, mesh->points[triangle->v[1]].xy));
    double len12 = sqrt(tm_distance_squared(mesh->points[triangle->v[1]].xy, mesh->points[triangle->v[2]].xy));
    double len20 = sqrt(tm_distance_squared(mesh->points[triangle->v[2]].xy, mesh->points[triangle->v[0]].xy));
    double shortest = len01;
    double area = fabs(orient2d(
        (REAL *) mesh->points[triangle->v[0]].xy,
        (REAL *) mesh->points[triangle->v[1]].xy,
        (REAL *) mesh->points[triangle->v[2]].xy
    )) * 0.5;
    double circumradius;
    double ratio;

    if (len12 < shortest) {
        shortest = len12;
    }
    if (len20 < shortest) {
        shortest = len20;
    }
    if (shortest <= 0.0 || area <= 0.0) {
        return 1;
    }

    circumradius = (len01 * len12 * len20) / (4.0 * area);
    ratio = circumradius / shortest;
    return ratio > beta * (1.0 + 1e-12);
}

TMStatus tm_validate_quality_mesh(const TMMesh *mesh, double min_angle_deg, TMValidationReport *out_report)
{
    TMStatus status;
    double beta;
    size_t segment_index;
    size_t tri_index;

    if (mesh == NULL || out_report == NULL || min_angle_deg <= 0.0 || min_angle_deg >= 90.0) {
        return TM_ERR_INTERNAL;
    }

    status = tm_validate_mesh(mesh, 1, out_report);
    beta = 1.0 / (2.0 * sin(min_angle_deg * tm_pi / 180.0));

    for (segment_index = 0; segment_index < mesh->segment_count; ++segment_index) {
        int encroaching_point = -1;

        if (!mesh->segments[segment_index].live) {
            continue;
        }
        if (mesh->segments[segment_index].is_protected) {
            continue;
        }
        if (tm_segment_is_encroached(mesh, segment_index, &encroaching_point)) {
            if (out_report->first_encroached_segment == (size_t) -1) {
                out_report->first_encroached_segment = segment_index;
                out_report->first_encroaching_point = encroaching_point;
            }
            out_report->encroached_segment_errors += 1;
        }
    }

    for (tri_index = 0; tri_index < mesh->triangle_count; ++tri_index) {
        if (tm_triangle_is_quality_exempt(mesh, (int) tri_index)) {
            out_report->exempt_triangle_count += 1;
            continue;
        }
        if (tm_triangle_is_bad(mesh, tri_index, beta)) {
            if (out_report->first_bad_triangle == (size_t) -1) {
                out_report->first_bad_triangle = tri_index;
            }
            out_report->bad_triangle_errors += 1;
        }
    }

    if (status != TM_OK ||
        out_report->encroached_segment_errors != 0 ||
        out_report->bad_triangle_errors != 0) {
        return TM_ERR_INVALID_MESH;
    }

    return TM_OK;
}
