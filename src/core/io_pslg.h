#ifndef DTCC_MESHER_IO_PSLG_H
#define DTCC_MESHER_IO_PSLG_H

#include <stddef.h>

#include "mesh.h"

typedef struct {
    TMPoint *points;
    size_t point_count;
    TMSegment *segments;
    size_t segment_count;
    double (*holes)[2];
    size_t hole_count;
} TMPSLG;

TMStatus tm_read_pslg_file(const char *path, TMPSLG *out_pslg);
TMStatus tm_validate_pslg(const TMPSLG *pslg);
TMStatus tm_validate_segment_graph(
    const TMPoint *points,
    size_t point_count,
    const TMSegment *segments,
    size_t segment_count
);
const char *tm_last_pslg_error_detail(void);
void tm_clear_pslg_error_detail(void);
void tm_set_pslg_error_detail(const char *format, ...);
void tm_free_pslg(TMPSLG *pslg);

#endif
