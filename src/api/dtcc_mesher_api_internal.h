#ifndef DTCC_MESHER_API_INTERNAL_H
#define DTCC_MESHER_API_INTERNAL_H

#include <stdarg.h>

#include "dtcc_mesher/dtcc_mesher.h"

#include "mesh.h"

dtcc_mesher_status dtcc_mesher_api_map_status(TMStatus status);
void dtcc_mesher_api_clear_error(dtcc_mesher_error *error);
void dtcc_mesher_api_set_error(dtcc_mesher_error *error, dtcc_mesher_status status, const char *format, ...);
void dtcc_mesher_api_set_error_va(dtcc_mesher_error *error, dtcc_mesher_status status, const char *format, va_list args);

#endif
