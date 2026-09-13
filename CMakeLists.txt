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
                   src/fusion_cross_sections.cpp src/fusion_kinetics.cpp
                   src/fusion_coulomb.cpp src/fusion_beam.cpp src/fusion_rate_model.cpp
                   src/fusion_products.cpp src/fusion_laboratory.cpp
                   src/fusion_nuclear_coulomb.cpp
                   src/fusion_alpha_amplitudes.cpp
                   src/fusion_alpha_spectrum.cpp src/fusion_pb_birth.cpp
                   src/fusion_target_burn.cpp
                   src/fusion_transfer_rate.cpp src/fusion_two_component.cpp
                   src/fusion_handoff.cpp src/fusion_thermal_burn.cpp src/fusion_nuclear_data.cpp src/fusion_source_state.cpp)
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

    add_library(fusion_kinetics_fortran STATIC
                fortran/fusion_kinetics_fortran.f90)
    add_library(pb11::fusion_kinetics_fortran ALIAS fusion_kinetics_fortran)
    set_target_properties(
        fusion_kinetics_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_kinetics_fortran
    )
    target_include_directories(
        fusion_kinetics_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_kinetics_fortran PUBLIC pb11::pb11)

    add_library(fusion_two_component_fortran STATIC
                fortran/fusion_two_component_fortran.f90)
    add_library(pb11::fusion_two_component_fortran ALIAS fusion_two_component_fortran)
    set_target_properties(fusion_two_component_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_two_component_fortran)
    target_include_directories(fusion_two_component_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_two_component_fortran PUBLIC
        pb11::pb11 fusion_kinetics_fortran)

    add_library(fusion_handoff_fortran STATIC fortran/fusion_handoff_fortran.f90)
    add_library(pb11::fusion_handoff_fortran ALIAS fusion_handoff_fortran)
    set_target_properties(fusion_handoff_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_handoff_fortran)
    target_include_directories(fusion_handoff_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_handoff_fortran PUBLIC pb11::pb11)

    add_library(fusion_pb_birth_fortran STATIC fortran/fusion_pb_birth_fortran.f90)
    add_library(pb11::fusion_pb_birth_fortran ALIAS fusion_pb_birth_fortran)
    set_target_properties(fusion_pb_birth_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_pb_birth_fortran)
    target_include_directories(fusion_pb_birth_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_pb_birth_fortran PUBLIC pb11::pb11)

    add_library(fusion_fsci_fortran STATIC fortran/fusion_fsci_fortran.f90)
    add_library(pb11::fusion_fsci_fortran ALIAS fusion_fsci_fortran)
    set_target_properties(fusion_fsci_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_fsci_fortran)
    target_include_directories(fusion_fsci_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_fsci_fortran PUBLIC pb11::pb11 fusion_spectrum_fortran)

    add_library(fusion_source_state_fortran STATIC fortran/fusion_source_state_fortran.f90)
    add_library(pb11::fusion_source_state_fortran ALIAS fusion_source_state_fortran)
    set_target_properties(fusion_source_state_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_source_state_fortran)
    target_include_directories(fusion_source_state_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_source_state_fortran PUBLIC pb11::pb11)

    add_library(fusion_nuclear_data_fortran STATIC fortran/fusion_nuclear_data_fortran.f90)
    add_library(pb11::fusion_nuclear_data_fortran ALIAS fusion_nuclear_data_fortran)
    set_target_properties(fusion_nuclear_data_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_nuclear_data_fortran)
    target_include_directories(fusion_nuclear_data_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_nuclear_data_fortran PUBLIC pb11::pb11 fusion_laboratory_fortran)

    add_library(fusion_thermal_burn_fortran STATIC fortran/fusion_thermal_burn_fortran.f90)
    add_library(pb11::fusion_thermal_burn_fortran ALIAS fusion_thermal_burn_fortran)
    set_target_properties(fusion_thermal_burn_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_thermal_burn_fortran)
    target_include_directories(fusion_thermal_burn_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_thermal_burn_fortran PUBLIC pb11::pb11)

    add_library(fusion_beam_fortran STATIC fortran/fusion_beam_fortran.f90)
    add_library(pb11::fusion_beam_fortran ALIAS fusion_beam_fortran)
    set_target_properties(
        fusion_beam_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_beam_fortran
    )
    target_include_directories(
        fusion_beam_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_beam_fortran PUBLIC pb11::pb11)

    add_library(fusion_rate_model_fortran STATIC fortran/fusion_rate_model_fortran.f90)
    add_library(pb11::fusion_rate_model_fortran ALIAS fusion_rate_model_fortran)
    set_target_properties(fusion_rate_model_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
        EXPORT_NAME fusion_rate_model_fortran)
    target_include_directories(fusion_rate_model_fortran INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>")
    target_link_libraries(fusion_rate_model_fortran PUBLIC pb11::pb11 fusion_beam_fortran)

    add_library(fusion_target_burn_fortran STATIC
                fortran/fusion_target_burn_fortran.f90)
    add_library(pb11::fusion_target_burn_fortran ALIAS
                fusion_target_burn_fortran)
    set_target_properties(
        fusion_target_burn_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_target_burn_fortran
    )
    target_include_directories(
        fusion_target_burn_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_target_burn_fortran PUBLIC
                          pb11::pb11 fusion_beam_fortran)

    add_library(fusion_products_fortran STATIC
                fortran/fusion_products_fortran.f90)
    add_library(pb11::fusion_products_fortran ALIAS fusion_products_fortran)
    set_target_properties(
        fusion_products_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_products_fortran
    )
    target_include_directories(
        fusion_products_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_products_fortran PUBLIC pb11::pb11)

    add_library(fusion_laboratory_fortran STATIC
                fortran/fusion_laboratory_fortran.f90)
    add_library(pb11::fusion_laboratory_fortran ALIAS fusion_laboratory_fortran)
    set_target_properties(
        fusion_laboratory_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_laboratory_fortran
    )
    target_include_directories(
        fusion_laboratory_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_laboratory_fortran PUBLIC pb11::pb11)

    add_library(fusion_spectrum_fortran STATIC
                fortran/fusion_spectrum_fortran.f90)
    add_library(pb11::fusion_spectrum_fortran ALIAS fusion_spectrum_fortran)
    set_target_properties(
        fusion_spectrum_fortran
        PROPERTIES
            Fortran_MODULE_DIRECTORY
                "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod"
            EXPORT_NAME fusion_spectrum_fortran
    )
    target_include_directories(
        fusion_spectrum_fortran
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-mod>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/fortran/modules>"
    )
    target_link_libraries(fusion_spectrum_fortran PUBLIC pb11::pb11)
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

    add_executable(test_fusion_kinetics tests/test_fusion_kinetics.cpp)
    target_link_libraries(test_fusion_kinetics PRIVATE pb11::pb11)
    add_test(NAME fusion_kinetics_tests COMMAND test_fusion_kinetics)

    add_executable(test_fusion_coulomb tests/test_fusion_coulomb.cpp)
    target_link_libraries(test_fusion_coulomb PRIVATE pb11::pb11)
    add_test(NAME fusion_coulomb_tests COMMAND test_fusion_coulomb)

    add_executable(test_fusion_slowing tests/test_fusion_slowing.cpp)
    target_link_libraries(test_fusion_slowing PRIVATE pb11::pb11)
    add_test(NAME fusion_slowing_tests COMMAND test_fusion_slowing)

    add_executable(test_fusion_rate_model tests/test_fusion_rate_model.cpp)
    target_link_libraries(test_fusion_rate_model PRIVATE pb11::pb11)
    add_test(NAME fusion_rate_model_tests COMMAND test_fusion_rate_model)

    add_executable(test_fusion_beam tests/test_fusion_beam.cpp)
    target_link_libraries(test_fusion_beam PRIVATE pb11::pb11)
    add_test(NAME fusion_beam_tests COMMAND test_fusion_beam)

    add_executable(test_fusion_transfer_rate tests/test_fusion_transfer_rate.cpp)
    target_link_libraries(test_fusion_transfer_rate PRIVATE pb11::pb11)
    add_test(NAME fusion_transfer_rate_tests COMMAND test_fusion_transfer_rate)
    add_executable(test_fusion_two_component tests/test_fusion_two_component.cpp)
    target_link_libraries(test_fusion_two_component PRIVATE pb11::pb11)
    add_test(NAME fusion_two_component_tests COMMAND test_fusion_two_component)

    add_executable(test_fusion_two_component_slowing tools/study_finite_temperature_split.cpp)
    target_compile_definitions(test_fusion_two_component_slowing PRIVATE FUSION_STUDY_PUBLIC_API=1)
    target_link_libraries(test_fusion_two_component_slowing PRIVATE pb11::pb11)
    add_test(NAME fusion_two_component_slowing_tests
        COMMAND test_fusion_two_component_slowing 1000 1000
            "${CMAKE_CURRENT_BINARY_DIR}/two-component-slowing.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/two-component-slowing.log")

    add_executable(test_fusion_handoff tests/test_fusion_handoff.cpp)
    target_link_libraries(test_fusion_handoff PRIVATE pb11::pb11)
    add_test(NAME fusion_handoff_tests COMMAND test_fusion_handoff)

    add_executable(test_fusion_pb_birth tests/test_fusion_pb_birth.cpp)
    target_link_libraries(test_fusion_pb_birth PRIVATE pb11::pb11)
    add_test(NAME fusion_pb_birth_tests COMMAND test_fusion_pb_birth)

    add_executable(test_fusion_fsci tests/test_fusion_fsci.cpp)
    target_link_libraries(test_fusion_fsci PRIVATE pb11::pb11)
    add_test(NAME fusion_fsci_tests COMMAND test_fusion_fsci)

    add_executable(test_fusion_source_state tests/test_fusion_source_state.cpp)
    target_link_libraries(test_fusion_source_state PRIVATE pb11::pb11)
    add_test(NAME fusion_source_state_tests COMMAND test_fusion_source_state)

    add_executable(test_fusion_nuclear_data tests/test_fusion_nuclear_data.cpp)
    target_link_libraries(test_fusion_nuclear_data PRIVATE pb11::pb11)
    add_test(NAME fusion_nuclear_data_tests COMMAND test_fusion_nuclear_data)

    add_executable(study_nuclear_windows tools/study_nuclear_windows.cpp)
    target_include_directories(study_nuclear_windows PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(study_nuclear_windows PRIVATE pb11::pb11)

    add_executable(validate_nuclear_birth tools/validate_nuclear_birth.cpp)
    target_link_libraries(validate_nuclear_birth PRIVATE pb11::pb11)
    add_test(NAME fusion_nuclear_birth_four_momentum_tests COMMAND validate_nuclear_birth)

    add_executable(test_fusion_thermal_burn tests/test_fusion_thermal_burn.cpp)
    target_link_libraries(test_fusion_thermal_burn PRIVATE pb11::pb11)
    add_test(NAME fusion_thermal_burn_tests COMMAND test_fusion_thermal_burn)

    add_executable(test_fusion_target_burn tests/test_fusion_target_burn.cpp)
    target_link_libraries(test_fusion_target_burn PRIVATE pb11::pb11)
    add_test(NAME fusion_target_burn_tests COMMAND test_fusion_target_burn)

    add_executable(test_fusion_products tests/test_fusion_products.cpp)
    target_link_libraries(test_fusion_products PRIVATE pb11::pb11)
    add_test(NAME fusion_products_tests COMMAND test_fusion_products)

    add_executable(test_fusion_laboratory tests/test_fusion_laboratory.cpp)
    target_link_libraries(test_fusion_laboratory PRIVATE pb11::pb11)
    add_test(NAME fusion_laboratory_tests COMMAND test_fusion_laboratory)

    add_executable(test_fusion_spectrum tests/test_fusion_spectrum.cpp)
    target_link_libraries(test_fusion_spectrum PRIVATE pb11::pb11)
    add_test(NAME fusion_spectrum_tests COMMAND test_fusion_spectrum)

    add_executable(test_fusion_alpha_grid tests/test_fusion_alpha_grid.cpp)
    target_link_libraries(test_fusion_alpha_grid PRIVATE pb11::pb11)
    add_test(NAME fusion_alpha_grid_tests COMMAND test_fusion_alpha_grid)

    add_executable(validate_alpha_grid tools/validate_alpha_grid.cpp)
    target_link_libraries(validate_alpha_grid PRIVATE pb11::pb11)

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

        add_executable(test_fusion_kinetics_fortran
                       tests/test_fusion_kinetics_fortran.f90)
        set_target_properties(test_fusion_kinetics_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_kinetics_fortran PRIVATE
                              pb11::fusion_kinetics_fortran)
        add_test(NAME fusion_kinetics_fortran_binding_tests
                 COMMAND test_fusion_kinetics_fortran)

        add_executable(test_fusion_rate_model_fortran tests/test_fusion_rate_model_fortran.f90)
        set_target_properties(test_fusion_rate_model_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_rate_model_fortran PRIVATE pb11::fusion_rate_model_fortran)
        add_test(NAME fusion_rate_model_fortran_binding_tests COMMAND test_fusion_rate_model_fortran)

        add_executable(test_fusion_beam_fortran
                       tests/test_fusion_beam_fortran.f90)
        set_target_properties(test_fusion_beam_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_beam_fortran PRIVATE
                              pb11::fusion_beam_fortran)
        add_test(NAME fusion_beam_fortran_binding_tests
                 COMMAND test_fusion_beam_fortran)

        add_executable(test_fusion_two_component_fortran
                       tests/test_fusion_two_component_fortran.f90)
        set_target_properties(test_fusion_two_component_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_two_component_fortran PRIVATE
            pb11::fusion_two_component_fortran)
        add_test(NAME fusion_two_component_fortran_binding_tests
                 COMMAND test_fusion_two_component_fortran)

        add_executable(test_fusion_handoff_fortran tests/test_fusion_handoff_fortran.f90)
        set_target_properties(test_fusion_handoff_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_handoff_fortran PRIVATE pb11::fusion_handoff_fortran)
        add_test(NAME fusion_handoff_fortran_binding_tests COMMAND test_fusion_handoff_fortran)

        add_executable(test_fusion_pb_birth_fortran tests/test_fusion_pb_birth_fortran.f90)
        set_target_properties(test_fusion_pb_birth_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_pb_birth_fortran PRIVATE pb11::fusion_pb_birth_fortran)
        add_test(NAME fusion_pb_birth_fortran_binding_tests COMMAND test_fusion_pb_birth_fortran)

        add_executable(test_fusion_fsci_fortran tests/test_fusion_fsci_fortran.f90)
        set_target_properties(test_fusion_fsci_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_fsci_fortran PRIVATE pb11::fusion_fsci_fortran)
        add_test(NAME fusion_fsci_fortran_binding_tests COMMAND test_fusion_fsci_fortran)

        add_executable(test_fusion_source_state_fortran tests/test_fusion_source_state_fortran.f90)
        set_target_properties(test_fusion_source_state_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_source_state_fortran PRIVATE pb11::fusion_source_state_fortran)
        add_test(NAME fusion_source_state_fortran_binding_tests COMMAND test_fusion_source_state_fortran)

        add_executable(test_fusion_nuclear_data_fortran tests/test_fusion_nuclear_data_fortran.f90)
        set_target_properties(test_fusion_nuclear_data_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_nuclear_data_fortran PRIVATE pb11::fusion_nuclear_data_fortran)
        add_test(NAME fusion_nuclear_data_fortran_binding_tests COMMAND test_fusion_nuclear_data_fortran)

        add_executable(test_fusion_thermal_burn_fortran tests/test_fusion_thermal_burn_fortran.f90)
        set_target_properties(test_fusion_thermal_burn_fortran PROPERTIES
            LINKER_LANGUAGE Fortran
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-fortran-mod")
        target_link_libraries(test_fusion_thermal_burn_fortran PRIVATE pb11::fusion_thermal_burn_fortran)
        add_test(NAME fusion_thermal_burn_fortran_binding_tests COMMAND test_fusion_thermal_burn_fortran)

        add_executable(test_fusion_target_burn_fortran
                       tests/test_fusion_target_burn_fortran.f90)
        set_target_properties(test_fusion_target_burn_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_target_burn_fortran PRIVATE
                              pb11::fusion_target_burn_fortran)
        add_test(NAME fusion_target_burn_fortran_binding_tests
                 COMMAND test_fusion_target_burn_fortran)

        add_executable(test_fusion_products_fortran
                       tests/test_fusion_products_fortran.f90)
        set_target_properties(test_fusion_products_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_products_fortran PRIVATE
                              pb11::fusion_products_fortran)
        add_test(NAME fusion_products_fortran_binding_tests
                 COMMAND test_fusion_products_fortran)

        add_executable(test_fusion_laboratory_fortran
                       tests/test_fusion_laboratory_fortran.f90)
        set_target_properties(test_fusion_laboratory_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_laboratory_fortran PRIVATE
                              pb11::fusion_laboratory_fortran)
        add_test(NAME fusion_laboratory_fortran_binding_tests
                 COMMAND test_fusion_laboratory_fortran)

        add_executable(test_fusion_spectrum_fortran
                       tests/test_fusion_spectrum_fortran.f90)
        set_target_properties(test_fusion_spectrum_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_spectrum_fortran PRIVATE
                              pb11::fusion_spectrum_fortran)
        add_test(NAME fusion_spectrum_fortran_binding_tests
                 COMMAND test_fusion_spectrum_fortran)

        add_executable(test_fusion_alpha_grid_fortran
                       tests/test_fusion_alpha_grid_fortran.f90)
        set_target_properties(test_fusion_alpha_grid_fortran PROPERTIES
                              POSITION_INDEPENDENT_CODE ON
                              LINKER_LANGUAGE Fortran)
        target_link_libraries(test_fusion_alpha_grid_fortran PRIVATE
                              pb11::fusion_spectrum_fortran)
        add_test(NAME fusion_alpha_grid_fortran_binding_tests
                 COMMAND test_fusion_alpha_grid_fortran)
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
        TARGETS pb11_fortran fusion_fortran fusion_kinetics_fortran
                fusion_two_component_fortran fusion_handoff_fortran fusion_thermal_burn_fortran fusion_nuclear_data_fortran fusion_source_state_fortran
                fusion_beam_fortran fusion_rate_model_fortran fusion_target_burn_fortran fusion_fsci_fortran fusion_pb_birth_fortran
                fusion_products_fortran
                fusion_laboratory_fortran fusion_spectrum_fortran
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
          include/fusion_cross_sections.h include/fusion_kinetics.h
          include/fusion_coulomb.h include/fusion_beam.h include/fusion_rate_model.h
          include/fusion_target_burn.h include/fusion_two_component.h include/fusion_handoff.h
          include/fusion_thermal_burn.h include/fusion_nuclear_data.h include/fusion_source_state.h
          include/fusion_products.h include/fusion_laboratory.h
          include/fusion_nuclear_coulomb.h include/fusion_alpha_amplitudes.h
          include/fusion_alpha_spectrum.h include/fusion_pb_birth.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(
    EXPORT pb11Targets
    NAMESPACE pb11::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pb11_reaction"
)
