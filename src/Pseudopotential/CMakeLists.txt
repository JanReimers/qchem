add_library(qcPseudopotential)
target_sources(qcPseudopotential PRIVATE
    Imp/GTH_Potentials.C
    )
target_sources(qcPseudopotential PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    LocalPotential.C
    SeparablePotential.C
    Integrals_Pseudo.C
    GTH_Potentials.C
    )
target_include_directories(qcPseudopotential PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
set_property(TARGET qcPseudopotential PROPERTY CXX_STANDARD 20)
target_compile_options(qcPseudopotential PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcPseudopotential PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS} ${PseudoDataPath})
target_link_libraries(qcPseudopotential PUBLIC qcMath qcCommon qcStructure)
