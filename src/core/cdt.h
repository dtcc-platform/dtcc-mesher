#ifndef DTCC_MESHER_CDT_H
#define DTCC_MESHER_CDT_H

#include "io_pslg.h"
#include "mesh.h"

typedef enum {
    TM_ACUTE_MODE_SIMPLE = 0,
    TM_ACUTE_MODE_SHELL
} TMAcuteProtectionMode;

typedef struct {
    int verbose;
    int refine;
    int protect_acute_corners;
    TMAcuteProtectionMode acute_mode;
    double min_angle_deg;
    double max_area;
    double max_edge_length;
    double protect_angle_deg;
    size_t max_refinement_steps;
    size_t max_protection_levels;
} TMBuildOptions;

typedef struct {
    double xy[2];
    int marker;
} TMRegion;

TMStatus tm_build_pslg_mesh(const TMPSLG *pslg, const TMBuildOptions *options, TMMesh *out_mesh);
TMStatus tm_build_coverage_mesh(
    const TMPSLG *pslg,
    const TMRegion *regions,
    size_t region_count,
    const TMBuildOptions *options,
    TMMesh *out_mesh,
    int **out_triangle_markers
);

#endif
