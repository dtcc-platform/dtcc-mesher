#ifndef DTCC_MESHER_PUBLIC_DTCC_MESHER_IO_H
#define DTCC_MESHER_PUBLIC_DTCC_MESHER_IO_H

#include "dtcc_mesher.h"

#ifdef __cplusplus
extern "C" {
#endif

TM_EXPORT tm_status tm_read_domain_file(const char *path, tm_domain *out_domain, tm_error *out_error);
TM_EXPORT void tm_domain_free(tm_domain *domain);

TM_EXPORT tm_status tm_write_triangles(const tm_mesh *mesh, const char *path, tm_error *out_error);
TM_EXPORT tm_status tm_write_svg(const tm_mesh *mesh, const char *path, tm_error *out_error);
TM_EXPORT tm_status tm_write_quality_csv(const tm_mesh *mesh, const char *path, tm_error *out_error);
TM_EXPORT tm_status tm_write_quality_summary(const tm_mesh *mesh, const char *path, tm_error *out_error);

#ifdef __cplusplus
}
#endif

#endif
