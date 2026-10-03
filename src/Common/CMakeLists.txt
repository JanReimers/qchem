# hand make std.pcm:
#   clang++-20 -std=c++23 -stdlib=libc++     -Wno-reserved-identifier -Wno-reserved-module-identifier     --precompile -o std.pcm /usr/lib/llvm-20/share/libc++/v1/std.cppm
# qcCommon: the non-math residue (pure utility/data).  All math lives in qcMath now (../Math).
add_library(qcCommon)
target_sources(qcCommon PRIVATE  Imp/PeriodicTable.C Imp/Reporting.C Imp/Parallel.C Imp/Diagnostics.C Imp/RunPolicy.C Imp/HDF5.C)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    PeriodicTable.C
    Iterators.C
    Outcome.C
    stl_io.C
    stream.C
    Strings.C
    Reporting.C
    Parallel.C
    Diagnostics.C
    RunPolicy.C
    HDF5.C
    )

set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcCommon PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcCommon PUBLIC ${CommonDataPath})

target_include_directories(qcCommon PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include   # PeriodicTable JSON
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include       # Strings (table formatting)
)
# openblas: PinBlasToOneThread (qchem.Parallel) calls openblas_set_num_threads.  PRIVATE -- the
# dependency is one call inside this library, not part of its interface.  Naming the DSO here is
# deliberate: if the system BLAS is ever swapped back to netlib, THIS link fails loudly instead of
# the pin silently disappearing (which is how it was lost once already, in the 2026-08 machine
# migration -- see the rationale on the declaration).
target_link_libraries(qcCommon PRIVATE openblas)
# qchem.HDF5 (the saved-SCF-state carrier, CK-1): the SERIAL HDF5 C library, named explicitly -- this box
# also carries the OpenMPI build (libhdf5_openmpi), which FindHDF5 can pick up instead, and an MPI-linked
# HDF5 inside a serial program is the singleton-init hang waiting to happen (CLAUDE.md, oracle codes).
find_path(HDF5_SERIAL_INCLUDE hdf5.h PATHS /usr/include/hdf5/serial NO_DEFAULT_PATH REQUIRED)
find_library(HDF5_SERIAL_LIB NAMES hdf5_serial REQUIRED)
target_include_directories(qcCommon PRIVATE ${HDF5_SERIAL_INCLUDE})
target_link_libraries(qcCommon PRIVATE ${HDF5_SERIAL_LIB})

add_subdirectory(tests)