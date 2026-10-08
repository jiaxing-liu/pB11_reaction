find_package(Threads REQUIRED)
add_library(pb11_diagnostics STATIC capture_registry.cpp capture_trial.cpp trial_replay.cpp trial_output_digest.cpp)
add_library(pb11::diagnostics ALIAS pb11_diagnostics)
set_target_properties(pb11_diagnostics PROPERTIES EXPORT_NAME diagnostics)
target_compile_features(pb11_diagnostics PUBLIC cxx_std_17)
target_include_directories(pb11_diagnostics
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pb11_diagnostics PUBLIC pb11::pb11 PRIVATE Threads::Threads)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(pb11_diagnostics PRIVATE -Wall -Wextra -Wpedantic)
endif()
install(TARGETS pb11_diagnostics EXPORT pb11Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
add_executable(fusion_trial_replay trial_replay_main.cpp)
target_link_libraries(fusion_trial_replay PRIVATE pb11::diagnostics)
install(TARGETS fusion_trial_replay RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/fusion_capture_registry.h include/fusion_capture_trial.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(PB11_FORTRAN_AVAILABLE)
  add_library(fusion_capture_fortran STATIC fortran/fusion_capture_fortran.f90)
  add_library(pb11::capture_fortran ALIAS fusion_capture_fortran)
  set_target_properties(fusion_capture_fortran PROPERTIES
    EXPORT_NAME capture_fortran
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod")
  target_include_directories(fusion_capture_fortran INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>)
  target_link_libraries(fusion_capture_fortran PUBLIC pb11::diagnostics)
  install(TARGETS fusion_capture_fortran EXPORT pb11Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod/"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/fortran/modules
    FILES_MATCHING PATTERN "*.mod")
endif()

if(BUILD_TESTING)
  add_executable(test_diagnostic_sha256 tests/test_sha256_stream.cpp)
  target_include_directories(test_diagnostic_sha256 PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_compile_features(test_diagnostic_sha256 PRIVATE cxx_std_17)
  add_test(NAME diagnostic_sha256_stream COMMAND test_diagnostic_sha256)
  add_executable(test_diagnostic_registry tests/test_capture_registry.cpp)
  target_include_directories(test_diagnostic_registry PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_diagnostic_registry PRIVATE pb11::diagnostics)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_definitions(test_diagnostic_registry PRIVATE PB11_DIAGNOSTIC_DESTROY_OBSERVATION)
    target_link_options(test_diagnostic_registry PRIVATE
      -Wl,--wrap=fusion_c_birth_table_destroy -Wl,--wrap=fusion_c_beam_birth_table_destroy)
  endif()
  if(PB11_FORTRAN_AVAILABLE)
    enable_language(C)
    add_executable(test_diagnostic_trial_view_fortran
      tests/test_trial_view_fortran.f90 tests/trial_view_probe.c)
    target_link_libraries(test_diagnostic_trial_view_fortran PRIVATE pb11::capture_fortran)
    add_test(NAME diagnostic_trial_view_fortran COMMAND test_diagnostic_trial_view_fortran)
  endif()
  if(UNIX)
    add_executable(test_diagnostic_trial_snapshot tests/test_trial_snapshot.cpp)
    target_include_directories(test_diagnostic_trial_snapshot PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests)
    target_link_libraries(test_diagnostic_trial_snapshot PRIVATE pb11::diagnostics)
    # The tests deliberately require fresh output directories. Repeated CTest
    # invocations get new directories without deleting earlier evidence.
    add_test(NAME diagnostic_trial_snapshot COMMAND ${CMAKE_COMMAND}
      -DTEST_PROGRAM=$<TARGET_FILE:test_diagnostic_trial_snapshot>
      -DOUTPUT_PARENT=${CMAKE_CURRENT_BINARY_DIR}/trial-test-output
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fresh_test.cmake)
    add_executable(test_diagnostic_trial_tables tests/test_trial_tables.cpp)
    target_include_directories(test_diagnostic_trial_tables PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests)
    target_link_libraries(test_diagnostic_trial_tables PRIVATE pb11::diagnostics)
  endif()
  set(PB11_DIAGNOSTIC_THERMAL_TEST_FILE "" CACHE FILEPATH "Qualified packed thermal table for optional lifecycle test")
  set(PB11_DIAGNOSTIC_BEAM_TEST_FILE "" CACHE FILEPATH "Qualified packed beam table for optional lifecycle test")
  if(PB11_DIAGNOSTIC_THERMAL_TEST_FILE)
    if(PB11_DIAGNOSTIC_BEAM_TEST_FILE)
      add_test(NAME diagnostic_registry_lifecycle COMMAND test_diagnostic_registry
        ${PB11_DIAGNOSTIC_THERMAL_TEST_FILE} ${PB11_DIAGNOSTIC_BEAM_TEST_FILE})
      if(UNIX)
        add_test(NAME diagnostic_trial_table_identity COMMAND ${CMAKE_COMMAND}
          -DTEST_PROGRAM=$<TARGET_FILE:test_diagnostic_trial_tables>
          -DOUTPUT_PARENT=${CMAKE_CURRENT_BINARY_DIR}/trial-table-test-output
          -DTHERMAL_FILE=${PB11_DIAGNOSTIC_THERMAL_TEST_FILE}
          -DBEAM_FILE=${PB11_DIAGNOSTIC_BEAM_TEST_FILE}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fresh_test.cmake)
      endif()
    else()
      add_test(NAME diagnostic_registry_lifecycle COMMAND test_diagnostic_registry
        ${PB11_DIAGNOSTIC_THERMAL_TEST_FILE})
    endif()
  endif()
endif()
