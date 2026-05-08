add_library(qcFactory) 
target_sources(qcFactory PRIVATE )
target_sources(qcFactory PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  Factory.C
)
set_property(TARGET qcFactory PROPERTY CXX_STANDARD 20)
target_link_libraries(qcFactory PUBLIC qcAtom_BS1 qcAtom_BS qcMolecule_BS) 
target_compile_definitions(qcFactory PUBLIC ${CXX_COMMON_DEFINES}
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)

