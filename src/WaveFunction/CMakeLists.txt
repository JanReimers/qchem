add_library(qcWaveFunction) 
target_sources(qcWaveFunction PRIVATE Factory.C   Polarized_WF.C Composite_WF.C  Irrep_WF.C  UnPolarized_WF.C EnergyLevel.C)
target_include_directories(qcWaveFunction PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/WaveFunction
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
    
    )
set_property(TARGET qcWaveFunction PROPERTY CXX_STANDARD 20)
target_link_libraries(qcWaveFunction PUBLIC qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon) 
