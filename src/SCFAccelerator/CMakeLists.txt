add_library(qcSCFAccelerator) 
target_sources(qcSCFAccelerator PRIVATE SCFAccelerator_DIIS.C  SCFAccelerator_Null.C Factory.C)
target_include_directories(qcSCFAccelerator PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/SCFAccelerator  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    )