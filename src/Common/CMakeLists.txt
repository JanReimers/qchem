# hand make std.pcm:
#   clang++-20 -std=c++23 -stdlib=libc++     -Wno-reserved-identifier -Wno-reserved-module-identifier     --precompile -o std.pcm /usr/lib/llvm-20/share/libc++/v1/std.cppm
# qcCommon: the non-math residue (pure utility/data).  All math lives in qcMath now (../Math).
add_library(qcCommon)
target_sources(qcCommon PRIVATE  Imp/PeriodicTable.C Imp/Reporting.C Imp/Parallel.C)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    PeriodicTable.C
    Iterators.C
    Outcome.C
    stl_io.C
    stream.C
    Strings.C
    Reporting.C
    Parallel.C
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

add_subdirectory(tests)