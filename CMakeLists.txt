cmake_minimum_required(VERSION 3.20)

project(
  dtcc_mesher
  VERSION 0.1.0
  DESCRIPTION "MIT Licensed Two-Dimensional Quality Mesh Generator"
  LANGUAGES C CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(DTCC_MESHER_BUILD_SHARED "Build DTCC Mesher as a shared library" OFF)
option(DTCC_MESHER_BUILD_CLI "Build the DTCC Mesher CLI" ON)
option(DTCC_MESHER_BUILD_PYTHON "Build Python bindings" OFF)
option(DTCC_MESHER_BUILD_TESTS "Build tests" ON)
option(DTCC_MESHER_BUILD_DEMOS "Build demos" ON)

set(DTCC_MESHER_LIBRARY_TYPE STATIC)
if(DTCC_MESHER_BUILD_SHARED)
  set(DTCC_MESHER_LIBRARY_TYPE SHARED)
endif()

set(DTCC_MESHER_CORE_SOURCES
  src/core/mesh.c
  src/core/io_pslg.c
  src/core/cdt.c
  src/core/report.c
  src/core/validate.c
  src/third_party/predicates.c
  src/api/dtcc_mesher.c
  src/api/dtcc_mesher_io.c
)

add_library(dtcc_mesher ${DTCC_MESHER_LIBRARY_TYPE} ${DTCC_MESHER_CORE_SOURCES})
add_library(dtcc_mesher::dtcc_mesher ALIAS dtcc_mesher)

target_compile_features(dtcc_mesher PUBLIC c_std_11)
target_include_directories(dtcc_mesher
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party
    ${CMAKE_CURRENT_SOURCE_DIR}/src/api
)

if(MSVC)
  target_compile_options(dtcc_mesher PRIVATE /W4)
else()
  target_compile_options(dtcc_mesher PRIVATE -Wall -Wextra -pedantic)
  set_source_files_properties(
    src/third_party/predicates.c
    PROPERTIES
      COMPILE_OPTIONS "-Wno-old-style-definition;-Wno-strict-prototypes;-Wno-declaration-after-statement;-Wno-unused-function;-Wno-missing-prototypes"
  )
endif()

set_target_properties(dtcc_mesher PROPERTIES
  OUTPUT_NAME dtcc_mesher
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(NOT MSVC)
  target_link_libraries(dtcc_mesher PUBLIC m)
endif()

if(DTCC_MESHER_BUILD_CLI)
  add_executable(dtcc_mesher_cli cli/main.c)
  target_compile_features(dtcc_mesher_cli PRIVATE c_std_11)
  target_link_libraries(dtcc_mesher_cli PRIVATE dtcc_mesher)
  target_include_directories(dtcc_mesher_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if(MSVC)
    target_compile_options(dtcc_mesher_cli PRIVATE /W4)
  else()
    target_compile_options(dtcc_mesher_cli PRIVATE -Wall -Wextra -pedantic)
  endif()
  set_target_properties(dtcc_mesher_cli PROPERTIES OUTPUT_NAME dtcc_mesher)
endif()

if(DTCC_MESHER_BUILD_DEMOS)
  add_executable(dtcc_mesher_demo_basic demos/c/basic_generate.c)
  target_compile_features(dtcc_mesher_demo_basic PRIVATE c_std_11)
  target_link_libraries(dtcc_mesher_demo_basic PRIVATE dtcc_mesher)
  target_include_directories(dtcc_mesher_demo_basic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

  add_executable(dtcc_mesher_demo_quality demos/c/quality_report.c)
  target_compile_features(dtcc_mesher_demo_quality PRIVATE c_std_11)
  target_link_libraries(dtcc_mesher_demo_quality PRIVATE dtcc_mesher)
  target_include_directories(dtcc_mesher_demo_quality PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

  if(NOT MSVC)
    target_compile_options(dtcc_mesher_demo_basic PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(dtcc_mesher_demo_quality PRIVATE -Wall -Wextra -pedantic)
  endif()
endif()

if(DTCC_MESHER_BUILD_TESTS)
  include(CTest)
  enable_testing()

  add_executable(dtcc_mesher_algorithmic_tests tests/test_mesh.c)
  target_compile_features(dtcc_mesher_algorithmic_tests PRIVATE c_std_11)
  target_link_libraries(dtcc_mesher_algorithmic_tests PRIVATE dtcc_mesher)
  target_include_directories(
    dtcc_mesher_algorithmic_tests
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/core
      ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party
  )
  if(MSVC)
    target_compile_options(dtcc_mesher_algorithmic_tests PRIVATE /W4)
  else()
    target_compile_options(dtcc_mesher_algorithmic_tests PRIVATE -Wall -Wextra -pedantic)
  endif()
  add_test(NAME c_algorithmic COMMAND dtcc_mesher_algorithmic_tests)
  set_tests_properties(c_algorithmic PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(dtcc_mesher_c_api_tests tests/c/test_api.c)
  target_compile_features(dtcc_mesher_c_api_tests PRIVATE c_std_11)
  target_link_libraries(dtcc_mesher_c_api_tests PRIVATE dtcc_mesher)
  target_include_directories(dtcc_mesher_c_api_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_compile_definitions(dtcc_mesher_c_api_tests PRIVATE DTCC_MESHER_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  if(MSVC)
    target_compile_options(dtcc_mesher_c_api_tests PRIVATE /W4)
  else()
    target_compile_options(dtcc_mesher_c_api_tests PRIVATE -Wall -Wextra -pedantic)
  endif()
  add_test(NAME c_api COMMAND dtcc_mesher_c_api_tests)
  set_tests_properties(c_api PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  if(DTCC_MESHER_BUILD_CLI)
    add_executable(dtcc_mesher_c_api_dump tests/c/api_dump.c)
    target_compile_features(dtcc_mesher_c_api_dump PRIVATE c_std_11)
    target_link_libraries(dtcc_mesher_c_api_dump PRIVATE dtcc_mesher)
    target_include_directories(dtcc_mesher_c_api_dump PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_compile_definitions(dtcc_mesher_c_api_dump PRIVATE DTCC_MESHER_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    if(MSVC)
      target_compile_options(dtcc_mesher_c_api_dump PRIVATE /W4)
    else()
      target_compile_options(dtcc_mesher_c_api_dump PRIVATE -Wall -Wextra -pedantic)
    endif()

    add_test(
      NAME cli_generate
      COMMAND
        ${CMAKE_COMMAND}
        -DCLI=$<TARGET_FILE:dtcc_mesher_cli>
        -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/cases/square_hole_domain.pslg
        -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}/tests/cli_generate
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/run_cli_test.cmake
    )
  endif()

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(DTCC_MESHER_BUILD_PYTHON AND Python3_Interpreter_FOUND)
    add_test(
      NAME python_pytest
      COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tests/python
    )
    if(WIN32)
      set(_tm_py_path_sep ";")
    else()
      set(_tm_py_path_sep ":")
    endif()
    set_tests_properties(
      python_pytest
      PROPERTIES
        ENVIRONMENT
          "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/python;DTCC_MESHER_TEST_CLI=$<TARGET_FILE:dtcc_mesher_cli>;DTCC_MESHER_TEST_C_API_DUMP=$<TARGET_FILE:dtcc_mesher_c_api_dump>;DTCC_MESHER_TEST_CASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/cases"
    )
  endif()
endif()

if(DTCC_MESHER_BUILD_PYTHON)
  find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)

  pybind11_add_module(dtcc_mesher_python MODULE python/bindings.cpp)
  target_compile_features(dtcc_mesher_python PRIVATE cxx_std_17)
  target_link_libraries(dtcc_mesher_python PRIVATE dtcc_mesher)
  target_include_directories(dtcc_mesher_python PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  set_target_properties(dtcc_mesher_python PROPERTIES
    OUTPUT_NAME _core
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python/dtcc_mesher
  )
  add_custom_command(
    TARGET dtcc_mesher_python
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/python/dtcc_mesher
      ${CMAKE_CURRENT_BINARY_DIR}/python/dtcc_mesher
  )
  if(NOT WIN32)
    set_target_properties(dtcc_mesher_python PROPERTIES BUILD_RPATH "$ORIGIN")
  endif()
endif()

install(
  TARGETS dtcc_mesher
  EXPORT dtcc_mesherTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(DTCC_MESHER_BUILD_CLI)
  install(TARGETS dtcc_mesher_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(DTCC_MESHER_BUILD_PYTHON)
  if(SKBUILD)
    set(DTCC_MESHER_PYTHON_INSTALL_DIR "dtcc_mesher")
  else()
    set(DTCC_MESHER_PYTHON_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/dtcc_mesher/python/dtcc_mesher")
  endif()

  install(TARGETS dtcc_mesher_python LIBRARY DESTINATION ${DTCC_MESHER_PYTHON_INSTALL_DIR})
endif()

install(DIRECTORY include/dtcc_mesher DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  FILES
    LICENSE
    README.md
    THIRD_PARTY.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/dtcc_mesher
)

if(DTCC_MESHER_BUILD_PYTHON)
  install(
    DIRECTORY python/dtcc_mesher/
    DESTINATION ${DTCC_MESHER_PYTHON_INSTALL_DIR}
    FILES_MATCHING
      PATTERN "*.py"
      PATTERN "py.typed"
  )
endif()

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesherConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dtcc_mesherConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesherConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dtcc_mesher
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dtcc_mesher.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesher.pc
  @ONLY
)

install(
  EXPORT dtcc_mesherTargets
  FILE dtcc_mesherTargets.cmake
  NAMESPACE dtcc_mesher::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dtcc_mesher
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesherConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesherConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dtcc_mesher
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/dtcc_mesher.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
