add_library(qcAtom_BS) 
target_sources(qcAtom_BS PRIVATE 
    Imp/Factory.C
    Imp/IE_HF.C
    Imp/IBS.C
    Imp/BS_Evaluator.C
    Internal/Imp/AngularIntegrals.C
    Internal/Imp/ExponentGrouper.C
    Internal/Imp/Exponential_Evaluator.C
    Internal/Imp/Wigner3j.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Factory.C
    IE.C
    IBS.C
    IBS_Evaluator.C
    BS_Evaluator.C
    Rk.C
    Internal/AngularIntegrals.C
    Internal/Wigner3j.C
    Internal/ExponentGrouper.C
    Internal/Exponential_Evaluator.C
    Internal/FillPower.C
    Internal/PascalTriangle.C
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
add_subdirectory(BSpline)
add_subdirectory(Gaussian)
add_subdirectory(Slater)
