add_library(qcAtom_BS) 
target_sources(qcAtom_BS PRIVATE 
    Imp/IEC.C
    Imp/Factory.C
    Imp/IE_HF.C
    Imp/BS.C
    Imp/IBS.C
    Internal/Imp/AngularIntegrals.C
    Internal/Imp/BFGrouper.C
    Internal/Imp/Wigner3j.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Factory.C
    IE.C
    IEC.C
    IBS.C
    BS.C
    IBS_Evaluator.C
    BS_Evaluator.C
    Rk.C
    Internal/AngularIntegrals.C
    Internal/BFGrouper.C
    Internal/Wigner3j.C
)
target_include_directories(qcAtom_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
target_compile_definitions(qcAtom_BS PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_compile_options(qcAtom_BS PRIVATE -Wall -Wno-unknown-pragmas -Wno-cpp 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
set_property(TARGET qcAtom_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcAtom_BS PUBLIC qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
add_subdirectory(l)
add_subdirectory(kappa)
add_subdirectory(radial)
