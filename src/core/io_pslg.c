#include "io_pslg.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "predicates.h"

static char *tm_next_data_line(FILE *stream, char *buffer, size_t buffer_size)
{
    while (fgets(buffer, (int) buffer_size, stream) != NULL) {
        char *cursor = buffer;

        while (*cursor != '\0' && isspace((unsigned char) *cursor)) {
            cursor += 1;
        }
        if (*cursor == '\0' || *cursor == '\n' || *cursor == '#') {
            continue;
        }
        return cursor;
    }

    return NULL;
}

static int tm_point_equals_xy(const TMPoint *point, const double xy[2])
{
    return point->xy[0] == xy[0] && point->xy[1] == xy[1];
}

static double tm_orient_xy(const double a[2], const double b[2], const double c[2])
{
    return orient2d((REAL *) a, (REAL *) b, (REAL *) c);
}

static int tm_point_on_closed_segment(const double a[2], const double b[2], const double p[2])
{
    if (tm_orient_xy(a, b, p) != 0.0) {
        return 0;
    }

    return p[0] >= ((a[0] < b[0]) ? a[0] : b[0]) &&
           p[0] <= ((a[0] > b[0]) ? a[0] : b[0]) &&
           p[1] >= ((a[1] < b[1]) ? a[1] : b[1]) &&
           p[1] <= ((a[1] > b[1]) ? a[1] : b[1]);
}

static int tm_segments_share_endpoint(const TMSegment *first, const TMSegment *second, const TMPoint *points)
{
    int i;
    int j;

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            if (tm_point_equals_xy(&points[first->v[i]], points[second->v[j]].xy)) {
                return 1;
            }
        }
    }

    return 0;
}

static int tm_segments_intersect_closed(const double a[2], const double b[2], const double c[2], const double d[2])
{
    double o1 = tm_orient_xy(a, b, c);
    double o2 = tm_orient_xy(a, b, d);
    double o3 = tm_orient_xy(c, d, a);
    double o4 = tm_orient_xy(c, d, b);

    if (((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0)) &&
        ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))) {
        return 1;
    }

    if (o1 == 0.0 && tm_point_on_closed_segment(a, b, c)) {
        return 1;
    }
    if (o2 == 0.0 && tm_point_on_closed_segment(a, b, d)) {
        return 1;
    }
    if (o3 == 0.0 && tm_point_on_closed_segment(c, d, a)) {
        return 1;
    }
    if (o4 == 0.0 && tm_point_on_closed_segment(c, d, b)) {
        return 1;
    }

    return 0;
}

static TMStatus tm_validate_segments(const TMPoint *points, size_t point_count, const TMSegment *segments, size_t segment_count)
{
    size_t *degree = NULL;
    size_t first_index;
    TMStatus status = TM_OK;

    degree = (size_t *) calloc(point_count, sizeof(*degree));
    if (degree == NULL) {
        return TM_ERR_ALLOC;
    }

    for (first_index = 0; first_index < segment_count; ++first_index) {
        const TMSegment *first = &segments[first_index];
        const double *a = points[first->v[0]].xy;
        const double *b = points[first->v[1]].xy;
        size_t second_index;

        if (first->v[0] < 0 || first->v[1] < 0 ||
            (size_t) first->v[0] >= point_count || (size_t) first->v[1] >= point_count) {
            status = TM_ERR_INVALID_PSLG;
            goto cleanup;
        }

        if (first->v[0] == first->v[1] || tm_point_equals_xy(&points[first->v[0]], points[first->v[1]].xy)) {
            status = TM_ERR_INVALID_PSLG;
            goto cleanup;
        }

        degree[first->v[0]] += 1;
        degree[first->v[1]] += 1;

        for (second_index = first_index + 1; second_index < segment_count; ++second_index) {
            const TMSegment *second = &segments[second_index];
            const double *c = points[second->v[0]].xy;
            const double *d = points[second->v[1]].xy;

            if ((first->v[0] == second->v[0] && first->v[1] == second->v[1]) ||
                (first->v[0] == second->v[1] && first->v[1] == second->v[0])) {
                status = TM_ERR_INVALID_PSLG;
                goto cleanup;
            }

            if (tm_segments_share_endpoint(first, second, points)) {
                int endpoint_ok = 1;

                if (second->v[0] != first->v[0] && second->v[0] != first->v[1] &&
                    tm_point_on_closed_segment(a, b, c)) {
                    endpoint_ok = 0;
                }
                if (second->v[1] != first->v[0] && second->v[1] != first->v[1] &&
                    tm_point_on_closed_segment(a, b, d)) {
                    endpoint_ok = 0;
                }
                if (first->v[0] != second->v[0] && first->v[0] != second->v[1] &&
                    tm_point_on_closed_segment(c, d, a)) {
                    endpoint_ok = 0;
                }
                if (first->v[1] != second->v[0] && first->v[1] != second->v[1] &&
                    tm_point_on_closed_segment(c, d, b)) {
                    endpoint_ok = 0;
                }

                if (!endpoint_ok) {
                    status = TM_ERR_INVALID_PSLG;
                    goto cleanup;
                }

                continue;
            }

            if (tm_segments_intersect_closed(a, b, c, d)) {
                status = TM_ERR_INVALID_PSLG;
                goto cleanup;
            }
        }
    }

    for (first_index = 0; first_index < point_count; ++first_index) {
        if (degree[first_index] != 0 && degree[first_index] != 2) {
            status = TM_ERR_INVALID_PSLG;
            goto cleanup;
        }
    }

cleanup:
    free(degree);
    return status;
}

TMStatus tm_read_pslg_file(const char *path, TMPSLG *out_pslg)
{
    FILE *stream;
    char line[4096];
    char keyword[64];
    size_t count = 0;
    size_t i;
    TMStatus status;

    if (path == NULL || out_pslg == NULL) {
        return TM_ERR_INTERNAL;
    }

    memset(out_pslg, 0, sizeof(*out_pslg));

    stream = fopen(path, "r");
    if (stream == NULL) {
        return TM_ERR_IO;
    }

    if (tm_next_data_line(stream, line, sizeof(line)) == NULL || sscanf(line, "%63s %zu", keyword, &count) != 2 ||
        strcmp(keyword, "vertices") != 0) {
        fclose(stream);
        return TM_ERR_PARSE;
    }

    out_pslg->points = (TMPoint *) calloc(count, sizeof(*out_pslg->points));
    if (out_pslg->points == NULL) {
        fclose(stream);
        return TM_ERR_ALLOC;
    }
    out_pslg->point_count = count;

    for (i = 0; i < count; ++i) {
        double x;
        double y;

        if (tm_next_data_line(stream, line, sizeof(line)) == NULL || sscanf(line, "%lf %lf", &x, &y) != 2) {
            fclose(stream);
            tm_free_pslg(out_pslg);
            return TM_ERR_PARSE;
        }

        out_pslg->points[i].xy[0] = x;
        out_pslg->points[i].xy[1] = y;
        out_pslg->points[i].original_index = (int) i;
        out_pslg->points[i].kind = TM_VERTEX_INPUT;
        out_pslg->points[i].incident_triangle = -1;
    }

    if (tm_next_data_line(stream, line, sizeof(line)) == NULL || sscanf(line, "%63s %zu", keyword, &count) != 2 ||
        strcmp(keyword, "segments") != 0) {
        fclose(stream);
        tm_free_pslg(out_pslg);
        return TM_ERR_PARSE;
    }

    out_pslg->segments = (TMSegment *) calloc(count, sizeof(*out_pslg->segments));
    if (out_pslg->segments == NULL) {
        fclose(stream);
        tm_free_pslg(out_pslg);
        return TM_ERR_ALLOC;
    }
    out_pslg->segment_count = count;

    for (i = 0; i < count; ++i) {
        int a;
        int b;

        if (tm_next_data_line(stream, line, sizeof(line)) == NULL || sscanf(line, "%d %d", &a, &b) != 2) {
            fclose(stream);
            tm_free_pslg(out_pslg);
            return TM_ERR_PARSE;
        }

        out_pslg->segments[i].v[0] = a;
        out_pslg->segments[i].v[1] = b;
        out_pslg->segments[i].original_index = (int) i;
        out_pslg->segments[i].live = 1;
    }

    if (tm_next_data_line(stream, line, sizeof(line)) != NULL) {
        if (sscanf(line, "%63s %zu", keyword, &count) != 2 || strcmp(keyword, "holes") != 0) {
            fclose(stream);
            tm_free_pslg(out_pslg);
            return TM_ERR_PARSE;
        }

        out_pslg->holes = (double (*)[2]) calloc(count, sizeof(*out_pslg->holes));
        if (out_pslg->holes == NULL && count != 0) {
            fclose(stream);
            tm_free_pslg(out_pslg);
            return TM_ERR_ALLOC;
        }
        out_pslg->hole_count = count;

        for (i = 0; i < count; ++i) {
            if (tm_next_data_line(stream, line, sizeof(line)) == NULL ||
                sscanf(line, "%lf %lf", &out_pslg->holes[i][0], &out_pslg->holes[i][1]) != 2) {
                fclose(stream);
                tm_free_pslg(out_pslg);
                return TM_ERR_PARSE;
            }
        }
    }

    fclose(stream);

    tm_initialize();
    status = tm_validate_segments(out_pslg->points, out_pslg->point_count, out_pslg->segments, out_pslg->segment_count);
    if (status != TM_OK) {
        tm_free_pslg(out_pslg);
        return status;
    }

    return TM_OK;
}

TMStatus tm_validate_pslg(const TMPSLG *pslg)
{
    if (pslg == NULL || (pslg->point_count != 0 && pslg->points == NULL)) {
        return TM_ERR_INTERNAL;
    }

    if (pslg->segment_count != 0 && pslg->segments == NULL) {
        return TM_ERR_INTERNAL;
    }

    if (pslg->hole_count != 0 && pslg->holes == NULL) {
        return TM_ERR_INTERNAL;
    }

    tm_initialize();
    return tm_validate_segments(pslg->points, pslg->point_count, pslg->segments, pslg->segment_count);
}

void tm_free_pslg(TMPSLG *pslg)
{
    if (pslg == NULL) {
        return;
    }

    free(pslg->points);
    free(pslg->segments);
    free(pslg->holes);
    pslg->points = NULL;
    pslg->segments = NULL;
    pslg->holes = NULL;
    pslg->point_count = 0;
    pslg->segment_count = 0;
    pslg->hole_count = 0;
}
