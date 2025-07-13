
add_library(qcCommon) 
target_sources(qcCommon PRIVATE   pmstream.C   stl_io.C  UniqueIDImp.C)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Constants.cpp Factorials.C PeriodicTable.C )
target_include_directories(qcCommon PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_options(qcCommon PRIVATE -Wall -Wno-unknown-pragmas -fmodules-ts
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
