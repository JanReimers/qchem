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

# Was a PRIVATE copy of the Debug/Release -O flags, which silently diverged from the common set: it
# missed -g -fno-omit-frame-pointer (so this lib alone was absent from perf call graphs) and, under
# C++20 modules, any GLOBAL ISA flag then fails to build outright -- a .pcm and its consumer must agree
# on target features (measured 2026-08-20: "compiled with the target feature '+avx2' but the current
target_compile_definitions(qcAtom_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})

set_property(TARGET qcAtom_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcAtom_BS PUBLIC qcBasisSet qcElConfig qcStructure lapack)
add_subdirectory(Evaluators)
add_subdirectory(tests)
