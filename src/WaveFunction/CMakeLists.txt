add_library(qcWaveFunction) 
target_sources(qcWaveFunction PRIVATE Imp/Factory.C   Internal/Imp/PolarizedWF.C Internal/Imp/CompositeWF.C  Internal/Imp/IrrepWF.C  Internal/Imp/UnPolarizedWF.C )
target_sources(qcWaveFunction PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    WaveFunction.C
    SCFWaveFunction.C
    Factory.C
    Types.C
    Internal/PolarizedWF.C 
    Internal/CompositeWF.C  
    Internal/IrrepWF.C   
    Internal/UnPolarizedWF.C
)
target_include_directories(qcWaveFunction PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    )
set_property(TARGET qcWaveFunction PROPERTY CXX_STANDARD 20)
target_link_libraries(qcWaveFunction PUBLIC qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcBasisSet qcSymmetry qcStructure qcCommon)
target_compile_options(qcWaveFunction PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcWaveFunction PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
