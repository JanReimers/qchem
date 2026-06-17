add_library(qcFactory) 
target_sources(qcFactory PRIVATE )
target_sources(qcFactory PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  Factory.C
)
set_property(TARGET qcFactory PROPERTY CXX_STANDARD 20)
target_link_libraries(qcFactory PUBLIC qcMolecule_BS qcAtom_BS)
target_compile_options(qcFactory PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcFactory PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})

