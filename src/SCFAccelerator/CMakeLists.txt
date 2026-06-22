add_library(qcSCFAccelerator) 
target_sources(qcSCFAccelerator PRIVATE Internal/Imp/SCFAcceleratorDIIS.C  Internal/Imp/SCFAcceleratorNull.C Internal/Imp/SCFAcceleratorGDM.C Internal/Imp/SCFAcceleratorLadder.C Imp/Factory.C)
target_sources(qcSCFAccelerator PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    SCFAccelerator.C
    Factory.C
    Internal/SCFAcceleratorDIIS.C
    Internal/SCFAcceleratorNull.C
    Internal/SCFAcceleratorGDM.C
    Internal/SCFAcceleratorLadder.C
)

target_include_directories(qcSCFAccelerator PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    )
 
set_property(TARGET qcSCFAccelerator PROPERTY CXX_STANDARD 20)
target_include_directories(qcSCFAccelerator PRIVATE 
      ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_compile_options(qcSCFAccelerator PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcSCFAccelerator PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcSCFAccelerator PUBLIC qcLASolver qcSymmetry qcStructure qcMesh qcCommon) 