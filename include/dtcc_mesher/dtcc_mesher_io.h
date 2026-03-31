#ifndef DTCC_MESHER_PUBLIC_DTCC_MESHER_IO_H
#define DTCC_MESHER_PUBLIC_DTCC_MESHER_IO_H

#include "dtcc_mesher.h"

#ifdef __cplusplus
extern "C" {
#endif

DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_read_domain_file(
    const char *path,
    dtcc_mesher_domain *out_domain,
    dtcc_mesher_error *out_error
);
DTCC_MESHER_EXPORT void dtcc_mesher_domain_free(dtcc_mesher_domain *domain);

DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_write_triangles(
    const dtcc_mesher_mesh *mesh,
    const char *path,
    dtcc_mesher_error *out_error
);
DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_write_svg(
    const dtcc_mesher_mesh *mesh,
    const char *path,
    dtcc_mesher_error *out_error
);
DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_write_quality_csv(
    const dtcc_mesher_mesh *mesh,
    const char *path,
    dtcc_mesher_error *out_error
);
DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_write_quality_summary(
    const dtcc_mesher_mesh *mesh,
    const char *path,
    dtcc_mesher_error *out_error
);

#ifdef __cplusplus
}
#endif

#endif
