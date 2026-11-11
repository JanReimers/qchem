add_library(qcPseudopotential)
target_sources(qcPseudopotential PRIVATE
    Imp/GTH_Potentials.C
    Imp/PseudopotentialSet.C
    )
target_sources(qcPseudopotential PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    LocalPotential.C
    SeparablePotential.C
    GTH_Potentials.C
    PseudopotentialSet.C
    )
target_include_directories(qcPseudopotential PRIVATE
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BLAZE_DIR}
)
set_property(TARGET qcPseudopotential PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcPseudopotential PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS} ${PseudoDataPath})
target_link_libraries(qcPseudopotential PUBLIC qcBasisSet qcMath qcCommon qcStructure)   # qcBasisSet: the Orbital_PP_IBS argument faces this library IMPLEMENTS (V1.2)

add_subdirectory(tests)
