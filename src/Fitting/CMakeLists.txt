add_library(qcFitting) 
target_sources(qcFitting PRIVATE ConstrainedFF.C  FittedFunction.C  IntegralConstrainedFF.C)
target_include_directories(qcFitting PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Fitting  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
