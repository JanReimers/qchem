# pybind/CMakeLists.txt -- Python bindings (qchem_py).  Built only when
# -DQCHEM_PYBIND=ON (the parent guards the add_subdirectory), so the normal
# build is completely unaffected.
#
# Requires PIC core libs: configure the whole project with
#   -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DQCHEM_PYBIND=ON
#   -DPython_EXECUTABLE=<your venv python>

find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)

include(FetchContent)
FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG        v2.4.0
)
FetchContent_MakeAvailable(nanobind)

# -- bridge: a C++20 MODULE UNIT that imports qchem and exposes the flat C ABI ---
add_library(qchem_bridge STATIC)
target_sources(qchem_bridge PUBLIC FILE_SET bridge_mod TYPE CXX_MODULES FILES qchem_bridge.cpp)
set_target_properties(qchem_bridge PROPERTIES
    CXX_STANDARD 20
    POSITION_INDEPENDENT_CODE ON)
target_compile_options    (qchem_bridge PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qchem_bridge PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_include_directories(qchem_bridge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})   # qcb_api.h
target_include_directories(qchem_bridge PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include)  # nlohmann/json.hpp (in the GMF)
# Linking the qchem aggregate pulls every qc* lib AND their CXX_MODULES file-sets.
target_link_libraries(qchem_bridge PUBLIC qchem)

# -- nanobind module: pure POD header, no qchem modules -------------------------
nanobind_add_module(qchem_py qchem_py.cpp)
set_target_properties(qchem_py PROPERTIES CXX_STANDARD 20)
target_link_libraries(qchem_py PRIVATE qchem_bridge)
