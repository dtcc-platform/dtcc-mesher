#ifndef DTCC_MESHER_API_INTERNAL_H
#define DTCC_MESHER_API_INTERNAL_H

#include <stdarg.h>

#include "dtcc_mesher/dtcc_mesher.h"

#include "mesh.h"

tm_status tm_api_map_status(TMStatus status);
void tm_api_clear_error(tm_error *error);
void tm_api_set_error(tm_error *error, tm_status status, const char *format, ...);
void tm_api_set_error_va(tm_error *error, tm_status status, const char *format, va_list args);

#endif
