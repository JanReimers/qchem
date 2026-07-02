add_library(qcAtom_BS) 
target_sources(qcAtom_BS PRIVATE 
    Imp/Factory.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Factory.C
    IrrepBasisSet.C
    BasisSet.C
)
target_include_directories(qcAtom_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(qcAtom_BS PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-O0>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O2>
)
target_compile_definitions(qcAtom_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})

set_property(TARGET qcAtom_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcAtom_BS PUBLIC qcBasisSet qcElConfig qcStructure lapack)
add_subdirectory(Evaluators)
add_subdirectory(tests)
