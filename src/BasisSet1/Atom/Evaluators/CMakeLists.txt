# add_library(qcAtom_BS1) 
target_sources(qcAtom_BS1 PRIVATE 
    Imp/Evaluator.C
    Internal/Imp/AngularIntegrals.C
    Internal/Imp/ExponentGrouper.C
    Internal/Imp/Exponential_Evaluator.C
    Internal/Imp/Wigner3j.C
)
target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Evaluator.C
    Concepts.C
    Internal/Rk.C
    Internal/AngularIntegrals.C
    Internal/Wigner3j.C
    Internal/ExponentGrouper.C
    Internal/Exponential_Evaluator.C
    Internal/FillPower.C
    Internal/PascalTriangle.C
)
# target_include_directories(qcAtom_BS1 PRIVATE 
#     ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
#     ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
#     ${PROJECT_SOURCE_DIR}/submodules/json/single_include
#     ${PROJECT_SOURCE_DIR}/submodules/blaze
#     ${PROJECT_SOURCE_DIR}/src
# )

add_subdirectory(BSpline)
add_subdirectory(Gaussian)
add_subdirectory(Slater)
