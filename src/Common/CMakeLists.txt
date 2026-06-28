# hand make std.pcm:
#   clang++-20 -std=c++23 -stdlib=libc++     -Wno-reserved-identifier -Wno-reserved-module-identifier     --precompile -o std.pcm /usr/lib/llvm-20/share/libc++/v1/std.cppm
# qcCommon: the non-math residue (pure utility/data).  All math lives in qcMath now (../Math).
add_library(qcCommon)
target_sources(qcCommon PRIVATE  Imp/PeriodicTable.C)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    PeriodicTable.C
    Iterators.C
    stl_io.C
    stream.C
    Strings.C
    )

set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_options(qcCommon PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcCommon PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcCommon PUBLIC ${CommonDataPath})

target_include_directories(qcCommon PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include   # PeriodicTable JSON
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include       # Strings (table formatting)
)
target_link_libraries(qcCommon PUBLIC)

add_subdirectory(tests)