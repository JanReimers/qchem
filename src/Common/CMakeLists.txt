# hand make std.pcm:
#   clang++-20 -std=c++23 -stdlib=libc++     -Wno-reserved-identifier -Wno-reserved-module-identifier     --precompile -o std.pcm /usr/lib/llvm-20/share/libc++/v1/std.cppm
add_library(qcCommon) 
target_sources(qcCommon PRIVATE  Imp/PeriodicTable.C Imp/OldPeriodicTable.C)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Constants.cpp 
    Factorials.C 
    PeriodicTable.C 
    UniqueID.C  
    UniqueIDImp.C 
    intpow.C 
    Iterators.C
    stl_io.C
    stream.C
    Types.C
    Blaze.C
    vector3d.C
    matrix3d.C
    Strings.C
    Math.C
    )

set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_options(qcCommon PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcCommon PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS} COMMON_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/Common/Data\")
target_include_directories(qcCommon PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)
target_link_libraries(qcCommon PUBLIC)

add_subdirectory(tests)