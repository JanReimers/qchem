
add_library(qc_Common) 
target_sources(qc_Common PRIVATE DFTDefines.C  PeriodicTable.C  pmstream.C  ScalarFunction.C  stl_io.C  UniqueIDImp.C  VectorFunction.C)
target_include_directories(qc_Common PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)

