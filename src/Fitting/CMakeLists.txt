add_library(qcFitting) 
target_sources(qcFitting PRIVATE Imp/ConstrainedFF.C  Imp/FittedFunction.C  Imp/IntegralConstrainedFF.C)
target_sources(qcFitting PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    FittedFunction.C
    FittedFunctionClient.C
    FittedFunctionImp.C
)
target_include_directories(qcFitting PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Fitting  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
set_property(TARGET qcFitting PROPERTY CXX_STANDARD 20)
target_link_libraries(qcFitting PUBLIC qcBasisSet qcSymmetry qcCluster qcMesh qcCommon) 

