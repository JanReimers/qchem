add_library(qcChargeDensity) 
target_sources(qcChargeDensity PRIVATE 
    Imp/ChargeDensity.C  
    Imp/CompositeCD.C  
    Imp/Factory.C
    Internal/Imp/IrrepCD.C
    Internal/Imp/FittedCDImp.C
    Internal/Imp/PolarizedCDImp.C
    Imp/Seed.C
    Imp/AtomicDensity.C
    Imp/CompositeFittedCD.C
    Imp/FourierSeedCD.C
)
target_sources(qcChargeDensity PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    ChargeDensity.C
    FittedCD.C
    FourierDensity.C
    CompositeCD.C
    Factory.C
    Seed.C
    AtomicDensity.C
    CompositeFittedCD.C
    FourierSeedCD.C
    Types.C
    Internal/IrrepCD.C
    Internal/FittedCDImp.C
    Internal/PolarizedCDImp.C
)

set_property(TARGET qcChargeDensity PROPERTY CXX_STANDARD 20)
target_include_directories(qcChargeDensity PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include)
target_link_libraries(qcChargeDensity PUBLIC qcFitting qcBasisSet qcSymmetry qcStructure qcMesh qcMath qcCommon)
target_compile_options(qcChargeDensity PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcChargeDensity PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcChargeDensity PRIVATE ${ChargeDensityDataPath})
