cmake_minimum_required(VERSION 3.16)

project(pb11_reaction VERSION 1.0.0 LANGUAGES CXX)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost 1.70 REQUIRED)

add_library(pb11 src/pb11.cpp)
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

include(CTest)
if(BUILD_TESTING)
    add_executable(test_pb11 tests/test_pb11.cpp)
    target_link_libraries(test_pb11 PRIVATE pb11::pb11)
    add_test(NAME pb11_unit_tests COMMAND test_pb11)
endif()
