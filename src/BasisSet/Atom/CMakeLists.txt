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
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/src
)

target_compile_definitions(qcAtom_BS PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_compile_options(qcAtom_BS PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>: -O2 -pg> #need to prevent -pg from appearing in fortran compile lines.
)

set_property(TARGET qcAtom_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcAtom_BS PUBLIC qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
add_subdirectory(Evaluators)
