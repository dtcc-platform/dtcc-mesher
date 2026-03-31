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
    int use_offcenters;
    int protect_acute_corners;
    TMAcuteProtectionMode acute_mode;
    double min_angle_deg;
    double protect_angle_deg;
    size_t max_refinement_steps;
    size_t max_protection_levels;
} TMBuildOptions;

TMStatus tm_build_pslg_mesh(const TMPSLG *pslg, const TMBuildOptions *options, TMMesh *out_mesh);

#endif
