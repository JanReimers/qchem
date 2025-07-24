# hand make std.pcm:
#   clang++-20 -std=c++23 -stdlib=libc++     -Wno-reserved-identifier -Wno-reserved-module-identifier     --precompile -o std.pcm /usr/lib/llvm-20/share/libc++/v1/std.cppm
add_library(qcCommon) 
target_sources(qcCommon PRIVATE  )
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
    )

target_include_directories(qcCommon PRIVATE )
set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_options(qcCommon PRIVATE -Wall -Wno-unknown-pragmas 
    $<$<CONFIG:Debug>:-g -O0 >
    $<$<CONFIG:Release>: -O2 >
)
target_link_libraries(qcCommon PUBLIC omlm )