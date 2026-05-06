add_library(qcChargeDensity) 
target_sources(qcChargeDensity PRIVATE 
    Imp/ChargeDensity.C  
    Imp/CompositeCD.C  
    Imp/Factory.C
    Internal/Imp/IrrepCD.C  
    Internal/Imp/FittedCDImp.C  
    Internal/Imp/PolarizedCDImp.C 
)
target_sources(qcChargeDensity PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    ChargeDensity.C 
    FittedCD.C 
    CompositeCD.C 
    Factory.C 
    Internal/IrrepCD.C 
    Internal/FittedCDImp.C 
    Internal/PolarizedCDImp.C
)

set_property(TARGET qcChargeDensity PROPERTY CXX_STANDARD 20)
target_include_directories(qcChargeDensity PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcChargeDensity PUBLIC qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon) 
