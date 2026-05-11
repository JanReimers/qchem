add_library(qcAtom_BS1) 
target_sources(qcAtom_BS1 PRIVATE 
    Imp/Factory.C
)
target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Factory.C
    IE.C
    IrrepBasisSet.C
    BasisSet.C
)
target_include_directories(qcAtom_BS1 PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/src
)

target_compile_definitions(qcAtom_BS1 PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_compile_options(qcAtom_BS1 PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>: -O2 -pg> #need to prevent -pg from appearing in fortran compile lines.
)

set_property(TARGET qcAtom_BS1 PROPERTY CXX_STANDARD 20)
target_link_libraries(qcAtom_BS1 PUBLIC qcBasisSet1 qcSymmetry qcCluster qcMesh qcCommon)
add_subdirectory(Evaluators)
