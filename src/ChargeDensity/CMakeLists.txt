add_library(qcChargeDensity) 
target_sources(qcChargeDensity PRIVATE Imp/ChargeDensity.C  Imp/CompositeCD.C  Imp/FittedCDImp.C  Imp/IrrepCDImp.C  Imp/PolarizedCDImp.C Imp/Factory.C)
target_sources(qcChargeDensity PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES ChargeDensity.C FittedCD.C CompositeCD.C Factory.C Imp/IrrepCD.C Imp/FittedCD.C Imp/PolarizedCD.C)

set_property(TARGET qcChargeDensity PROPERTY CXX_STANDARD 20)
target_link_libraries(qcChargeDensity PUBLIC qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon) 
