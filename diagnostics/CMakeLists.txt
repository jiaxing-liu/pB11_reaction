find_package(Threads REQUIRED)
add_library(pb11_diagnostics STATIC capture_registry.cpp)
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
install(FILES include/fusion_capture_registry.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
  set(PB11_DIAGNOSTIC_THERMAL_TEST_FILE "" CACHE FILEPATH "Qualified packed thermal table for optional lifecycle test")
  set(PB11_DIAGNOSTIC_BEAM_TEST_FILE "" CACHE FILEPATH "Qualified packed beam table for optional lifecycle test")
  if(PB11_DIAGNOSTIC_THERMAL_TEST_FILE)
    if(PB11_DIAGNOSTIC_BEAM_TEST_FILE)
      add_test(NAME diagnostic_registry_lifecycle COMMAND test_diagnostic_registry
        ${PB11_DIAGNOSTIC_THERMAL_TEST_FILE} ${PB11_DIAGNOSTIC_BEAM_TEST_FILE})
    else()
      add_test(NAME diagnostic_registry_lifecycle COMMAND test_diagnostic_registry
        ${PB11_DIAGNOSTIC_THERMAL_TEST_FILE})
    endif()
  endif()
endif()
