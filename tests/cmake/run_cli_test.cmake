file(MAKE_DIRECTORY "${OUTDIR}")

execute_process(
  COMMAND "${CLI}" "${INPUT}" "${OUTDIR}/mesh"
  RESULT_VARIABLE cli_result
)

if(NOT cli_result EQUAL 0)
  message(FATAL_ERROR "dtcc_mesher CLI test failed with exit code ${cli_result}")
endif()

foreach(path
  "${OUTDIR}/mesh.tri"
  "${OUTDIR}/mesh.svg"
  "${OUTDIR}/mesh.metrics.csv"
  "${OUTDIR}/mesh.summary.txt")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "missing CLI output ${path}")
  endif()
endforeach()

file(READ "${OUTDIR}/mesh.summary.txt" summary_text)

if(summary_text MATCHES "triangle_count=")
else()
  message(FATAL_ERROR "summary file does not contain triangle_count")
endif()
