cmake_minimum_required(VERSION 3.16)

project(pb11_reaction VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

option(PB11_BUILD_FORTRAN
       "Build the ISO_C_BINDING Fortran module when a compiler is available"
       ON)

set(PB11_FORTRAN_AVAILABLE OFF)
if(PB11_BUILD_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        set(PB11_FORTRAN_AVAILABLE ON)
    else()
        message(WARNING
                "Fortran compiler not found; pb11_fortran target and its test "
                "will be disabled")
    endif()
endif()

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost 1.70 REQUIRED)

add_library(pb11 src/pb11.cpp src/pb11_c.cpp src/pb11_source.cpp
                   src/fusion_network.cpp src/fusion_rates.cpp
                   src/fusion_cross_sections.cpp)
add_library(pb11::pb11 ALIAS pb11)

target_compile_features(pb11 PUBLIC cxx_std_17)
target_include_directories(
    pb11
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${Boost_INCLUDE_DIRS}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pb11 PRIVATE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(pb11 PRIVATE /W4)
endif()

add_executable(benchmark_pb11 examples/benchmark_pb11.cpp)
target_link_libraries(benchmark_pb11 PRIVATE pb11::pb11)
add_executable(thermal_network examples/thermal_network.cpp)
target_link_libraries(thermal_network PRIVATE pb11::pb11)

if(PB11_FORTRAN_AVAILABLE)
    add_library(pb11_fortran STATIC fortran/pb11_fortran.f90)
    add_library(pb11::fortran ALIAS pb11_fortran)
    set_target_properties(
        pb11_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fortran
    )
    target_include_directories(
        pb11_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(pb11_fortran PUBLIC pb11::pb11)

    add_library(fusion_fortran STATIC fortran/fusion_fortran.f90)
    add_library(pb11::fusion_fortran ALIAS fusion_fortran)
    set_target_properties(
        fusion_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_fortran
    )
    target_include_directories(
        fusion_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_fortran PUBLIC pb11::pb11)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(test_pb11 tests/test_pb11.cpp)
    target_link_libraries(test_pb11 PRIVATE pb11::pb11)
    add_test(NAME pb11_unit_tests COMMAND test_pb11)

    add_executable(test_pb11_c tests/test_pb11_c.cpp)
    target_link_libraries(test_pb11_c PRIVATE pb11::pb11)
    add_test(NAME pb11_c_abi_tests COMMAND test_pb11_c)

    add_executable(test_pb11_source tests/test_pb11_source.cpp)
    target_link_libraries(test_pb11_source PRIVATE pb11::pb11)
    add_test(NAME pb11_instant_source_tests COMMAND test_pb11_source)

    add_executable(test_fusion_network tests/test_fusion_network.cpp)
    target_link_libraries(test_fusion_network PRIVATE pb11::pb11)
    add_test(NAME fusion_network_tests COMMAND test_fusion_network)

    add_executable(test_fusion_rates tests/test_fusion_rates.cpp)
    target_link_libraries(test_fusion_rates PRIVATE pb11::pb11)
    add_test(NAME fusion_rates_tests COMMAND test_fusion_rates)

    add_executable(test_fusion_cross_sections
                   tests/test_fusion_cross_sections.cpp)
    target_link_libraries(test_fusion_cross_sections PRIVATE pb11::pb11)
    target_include_directories(test_fusion_cross_sections PRIVATE ${Boost_INCLUDE_DIRS})
    add_test(NAME fusion_cross_sections_tests COMMAND test_fusion_cross_sections)

    add_executable(test_fusion_closed_box tests/test_fusion_closed_box.cpp)
    target_link_libraries(test_fusion_closed_box PRIVATE pb11::pb11)
    add_test(NAME fusion_closed_box_tests COMMAND test_fusion_closed_box)

    if(PB11_FORTRAN_AVAILABLE)
        add_executable(test_pb11_fortran tests/test_pb11_fortran.f90)
        set_target_properties(test_pb11_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_pb11_fortran PRIVATE pb11::fortran)
        add_test(NAME pb11_fortran_binding_tests COMMAND test_pb11_fortran)

        add_executable(test_fusion_fortran tests/test_fusion_fortran.f90)
        set_target_properties(test_fusion_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_fortran PRIVATE pb11::fusion_fortran)
        add_test(NAME fusion_fortran_binding_tests COMMAND test_fusion_fortran)
    endif()
endif()

install(
    TARGETS pb11
    EXPORT pb11Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
if(PB11_FORTRAN_AVAILABLE)
    install(
        TARGETS pb11_fortran fusion_fortran
        EXPORT pb11Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    install(
        DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod/"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/fortran/modules"
        FILES_MATCHING PATTERN "*.mod"
    )
endif()
install(
    FILES include/pb11.hpp include/pb11_c.h include/pb11_source.h
          include/fusion_network.h include/fusion_rates.h
          include/fusion_cross_sections.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(
    EXPORT pb11Targets
    NAMESPACE pb11::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pb11_reaction"
)
