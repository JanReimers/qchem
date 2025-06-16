add_library(qcChargeDensity) 
target_sources(qcChargeDensity PRIVATE ChargeDensity.C  CompositeCD.C  FittedCD.C  IrrepCD.C  PolarizedCD.C)
target_include_directories(qcChargeDensity PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/ChargeDensity  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )