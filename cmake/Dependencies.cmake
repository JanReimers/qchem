# cmake/Dependencies.cmake -- THE PIN TABLE: every source dependency, at one exact commit, in one place
# (doc/SharedDependenciesPlan.md, ruled 2026-11-10; step 1).
#
# Each dependency is fetched as the GitHub archive of its pinned commit, verified by SHA-256 (the archive's content, so a
# re-tagged or rewritten upstream cannot change what builds), and placed under the build tree's _deps/.  The pins are the
# commits the retiring submodules held (each archive measured identical to its submodule checkout, 2026-11-10).  The FORKS
# (blaze, libxc) are the user's: change the fork, push, bump the pin + hash here.  To build against a LOCAL checkout while a
# fork patch is unpushed: -DFETCHCONTENT_SOURCE_DIR_<NAME>=/path (CMake's own override; NAME upper-case, e.g. LIBXC).
#
# Step 2 adds the shared install prefix (~/Code/deps/<compiler>/<dep>-<pin>, found first); this file stays the fallback for a
# tree with no prefix -- an outsider's clone, CI.
#
# What the project reads, after include():
#   QCHEM_BLAZE_DIR         blaze sources (headers; its cmake/ holds Blaze_Import)
#   QCHEM_JSON_INCLUDE      nlohmann/json single_include
#   QCHEM_BSPLINE_INCLUDE   BSplinebasis include
#   QCHEM_TABULATE_INCLUDE  tabulate include
#   QCHEM_POCKETFFT_DIR     pocketfft (pocketfft_hdronly.h)
#   QCHEM_GTEST_INCLUDE     googletest include; targets gtest / gtest_main (built here: a PIN, never the distro's -- ruling 2)
#   QCHEM_LIBXC_BUILD       libxc's BUILD dir (its generated headers); target xc
#   QCHEM_LIBCINT_SOURCE / QCHEM_LIBCINT_BUILD   libcint source and build dirs (include/ in each); target cint
include(FetchContent)

#  one row per dependency: name | repository (no .git) | commit | SHA-256 of the GitHub archive of that commit
set(QCHEM_PINS
    "blaze|https://github.com/JanReimers/blaze|7478a70b992dffb9a366a0be405730d1368f3eb9|4139b5d986ab724830912b95c66941e878116a7275bee16fb7583cde62b0b466"
    "json|https://github.com/nlohmann/json|568b708fd46deeb23a959381393a7564f1586588|283d800d754741c972e82d1a3bdccb10ef67904f090beed72bd6214da21e02c4"
    "BSplinebasis|https://github.com/okruz/BSplinebasis|678b28d226f781d4405ec20027489b520688d9a5|b599ba908ae068861d160ca59d929ac2aea2186132cdd3041d5f6df966d27c72"
    "tabulate|https://github.com/p-ranav/tabulate|b9026bf1e31dd561680e20cfda2a37716161aac1|edf6926366627047457e595d29c1c24b6fd44dfbda4d9b1df495254fd1f821f9"
    "pocketfft|https://github.com/mreineck/pocketfft|c90e55b3d529f8efa40ed01a20de22405f45fc65|2d90016da31a8cc1dd0ac0c89c1f2fb55df8f4ccc3f0386dcfbad09ac6d1c3ba"
    "googletest|https://github.com/google/googletest|175c1b55cfb3dbed519b94a370c083aac605009f|a21d61e995e222faa19f49e03cd53ca6985c44271f3fbdca2c56a2cd15e85543"
    "libxc|https://github.com/JanReimers/libxc|d00e5d0e2696cb15d7fcf6cf01dc447334cc70e2|7c1fc5ac3647acaf1bd13a844acd872e8079d52c7d4a2a295595a04059c18c78"
    "libcint|https://github.com/sunqm/libcint|c72d30662ba08048940c9000cc352a220c867712|c47663d737629742006f76895150aa9d5367ea4604adef6517f353b75075e0f3"
)

#! Declare every pin.  HEADER_ONLY ones are populated without add_subdirectory (SOURCE_SUBDIR points nowhere): their own
#! CMakeLists (tests, install rules, options) are not ours to configure.
set(_qchemHeaderOnly blaze json BSplinebasis tabulate pocketfft)
foreach(row IN LISTS QCHEM_PINS)
    string(REPLACE "|" ";" pin "${row}")   # one row: name | repository | commit | archive SHA-256
    list(GET pin 0 name)
    list(GET pin 1 repo)
    list(GET pin 2 commit)
    list(GET pin 3 sha)
    set(_sub "")
    if(name IN_LIST _qchemHeaderOnly)
        set(_sub SOURCE_SUBDIR qchem-headers-only-no-cmake)
    endif()
    FetchContent_Declare(${name}
        URL      ${repo}/archive/${commit}.tar.gz
        URL_HASH SHA256=${sha}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        ${_sub})
endforeach()

# The header-only sources.
FetchContent_MakeAvailable(blaze json BSplinebasis tabulate pocketfft)
set(QCHEM_BLAZE_DIR        ${blaze_SOURCE_DIR})
set(QCHEM_JSON_INCLUDE     ${json_SOURCE_DIR}/single_include)
set(QCHEM_BSPLINE_INCLUDE  ${bsplinebasis_SOURCE_DIR}/include)
set(QCHEM_TABULATE_INCLUDE ${tabulate_SOURCE_DIR}/include)
set(QCHEM_POCKETFFT_DIR    ${pocketfft_SOURCE_DIR})

# googletest: built here, at the pin (ruling 2: test discovery broke once on a googletest/CMake change -- never the distro's).
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock for this project" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "googletest: no install rules" FORCE)
FetchContent_MakeAvailable(googletest)
set(QCHEM_GTEST_INCLUDE ${googletest_SOURCE_DIR}/googletest/include)

# libxc / libcint are BUILT later (qchem_build_c_dependencies, after the project's targets), but src/ names their include
# directories before that: FetchContent's own layout says where they will be, and the macro checks it once they exist.
set(QCHEM_LIBXC_BUILD    ${FETCHCONTENT_BASE_DIR}/libxc-build)
set(QCHEM_LIBCINT_SOURCE ${FETCHCONTENT_BASE_DIR}/libcint-src)
set(QCHEM_LIBCINT_BUILD  ${FETCHCONTENT_BASE_DIR}/libcint-build)

#! The C libraries the project builds and links (libxc, libcint), made available where the root CMakeLists asks -- AFTER its
#! own targets, with BUILD_TESTING off (their tests are not ours), exactly as the retired add_subdirectory(submodules/...) did.
#! Both compile at plain -O3 (the root applies -march=native to C++ TUs only), which step 2's portable prefix relies on.
macro(qchem_build_c_dependencies)
    set(BUILD_TESTING OFF CACHE BOOL "Disable embedded dependency tests" FORCE)
    FetchContent_MakeAvailable(libxc)
    # libcint (Apache-2.0): static, no Fortran interface (we call the C API directly).
    set(ENABLE_STATIC ON  CACHE BOOL "Build libcint as a static library" FORCE)
    set(WITH_FORTRAN  OFF CACHE BOOL "libcint Fortran interface (unused)" FORCE)
    FetchContent_MakeAvailable(libcint)
    set_property(TARGET cint PROPERTY POSITION_INDEPENDENT_CODE ON)   # qcGaussian_BS is shared, so the static cint must be PIC
    set(BUILD_TESTING ON CACHE BOOL "Enable testing for this project" FORCE)
    foreach(_p "QCHEM_LIBXC_BUILD;libxc_BINARY_DIR" "QCHEM_LIBCINT_SOURCE;libcint_SOURCE_DIR" "QCHEM_LIBCINT_BUILD;libcint_BINARY_DIR")
        list(GET _p 0 _named)
        list(GET _p 1 _actual)
        if(NOT "${${_named}}" STREQUAL "${${_actual}}")
            message(FATAL_ERROR "Dependencies.cmake: ${_named} was announced as '${${_named}}' but FetchContent put it at '${${_actual}}'")
        endif()
    endforeach()
endmacro()
