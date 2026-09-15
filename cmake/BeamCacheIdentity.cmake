# Configure-time identity for code/data cached by the beam-table path.
#
# Call this after find_package(Boost REQUIRED), for example:
#
#   include("${CMAKE_CURRENT_LIST_DIR}/BeamCacheIdentity.cmake")
#   pb11_add_beam_cache_kernel_identity(pb11 "${PROJECT_SOURCE_DIR}")
#
# The library_root argument must contain src/ and include/.  The resulting
# header is generated below the caller's binary directory, so no source-tree
# or build-tree path becomes part of the identity value.

include_guard(GLOBAL)

set(_PB11_BEAM_CACHE_IDENTITY_TEMPLATE
    "${CMAKE_CURRENT_LIST_DIR}/fusion_cache_identity_internal.h.in")

function(pb11_add_beam_cache_kernel_identity target library_root)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR
      "pb11_add_beam_cache_kernel_identity: target does not exist: ${target}")
  endif()

  if(NOT DEFINED Boost_VERSION OR "${Boost_VERSION}" STREQUAL "")
    message(FATAL_ERROR
      "pb11_add_beam_cache_kernel_identity: call after find_package(Boost)")
  endif()

  get_filename_component(_root "${library_root}" ABSOLUTE
                         BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  if(NOT IS_DIRECTORY "${_root}")
    message(FATAL_ERROR
      "pb11_add_beam_cache_kernel_identity: library root is not a directory: ${library_root}")
  endif()

  # GLOB_RECURSE keeps the three required relative roots explicit while
  # CONFIGURE_DEPENDS makes additions/removals trigger a reconfigure.
  file(GLOB_RECURSE _relative_files
       LIST_DIRECTORIES false
       CONFIGURE_DEPENDS
       RELATIVE "${_root}"
       "${_root}/src/*.cpp"
       "${_root}/src/*.h"
       "${_root}/include/*.h"
       "${_root}/src/*.hpp" "${_root}/include/*.hpp"
       "${_root}/src/*.inc" "${_root}/include/*.inc")
  list(SORT _relative_files)

  if(NOT _relative_files)
    message(FATAL_ERROR
      "pb11_add_beam_cache_kernel_identity: no src/*.cpp, src/*.h, or include/*.h files under ${library_root}")
  endif()

  # The hash input is deliberately independent of absolute paths, build
  # paths, compiler flags, git state, and binary provenance.
  set(_hash_input "boost-version=${Boost_VERSION}\n")
  foreach(_relative_file IN LISTS _relative_files)
    string(REPLACE "\\" "/" _relative_file "${_relative_file}")
    set(_absolute_file "${_root}/${_relative_file}")
    file(SHA256 "${_absolute_file}" _file_sha256)

    # CONFIGURE_DEPENDS on the GLOB covers directory entries.  This property
    # covers edits to each currently matched file before the next build.
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 "${_absolute_file}")

    string(APPEND _hash_input
           "file=${_relative_file}\n"
           "sha256=${_file_sha256}\n")
  endforeach()
  string(SHA256 _kernel_identity "${_hash_input}")

  set(_generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
  file(MAKE_DIRECTORY "${_generated_dir}")
  set(PB11_BEAM_CACHE_KERNEL_IDENTITY "${_kernel_identity}")
  configure_file(
    "${_PB11_BEAM_CACHE_IDENTITY_TEMPLATE}"
    "${_generated_dir}/fusion_cache_identity_internal.h"
    @ONLY)

  target_include_directories("${target}" PRIVATE "${_generated_dir}")
endfunction()
