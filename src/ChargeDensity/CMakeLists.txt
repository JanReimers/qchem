add_library(qcChargeDensity) 
target_sources(qcChargeDensity PRIVATE Imp/ChargeDensity.C  Imp/CompositeCD.C  Imp/FittedCD.C  Imp/IrrepCDImp.C  PolarizedCD.C Factory.C)
target_sources(qcChargeDensity PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES ChargeDensity.C FittedCD.C CompositeCD.C Imp/IrrepCD.C)
target_include_directories(qcChargeDensity PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/ChargeDensity  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )

set_property(TARGET qcChargeDensity PROPERTY CXX_STANDARD 20)
target_link_libraries(qcChargeDensity PUBLIC qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon) 
