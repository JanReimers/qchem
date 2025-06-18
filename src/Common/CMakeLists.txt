
add_library(qcCommon) 
target_sources(qcCommon PRIVATE DFTDefines.C  Factorials.C PeriodicTable.C  pmstream.C   stl_io.C  UniqueIDImp.C )
target_include_directories(qcCommon PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
set_property(TARGET qcCommon PROPERTY CXX_STANDARD 20)
