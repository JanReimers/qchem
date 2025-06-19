
add_library(qcCommon) 
target_sources(qcCommon PRIVATE Factorials.C PeriodicTable.C  pmstream.C   stl_io.C  UniqueIDImp.C Constants.cpp)
target_sources(qcCommon PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Constants.cpp)
target_include_directories(qcCommon PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
target_compile_options(qcCommon PRIVATE -Wall -Wno-unknown-pragmas -fmodules-ts
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
