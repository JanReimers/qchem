# qcMath: the math leaf library (pulled out of qcCommon).  Independent of qcCommon.
#   - qchem.CMath  : <cmath> + Constants/IntPow/Factorials (the leaf the math core imports)
#   - qchem.Math   : lean umbrella for outside consumers (re-exports CMath)
#   - Types/Blaze, fixed Vector/Matrix 2D+3D, SpecialFunctions, FFT, FourierMap,
#     ScalarFunction, VectorFunction
add_library(qcMath)
target_sources(qcMath PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Constants.cpp
    intpow.C
    Factorials.C
    CMath.C
    Math.C
    Blaze.C
    Types.C
    vector2d.C
    vector3d.C
    matrix2d.C
    matrix3d.C
    SpecialFunctions.C
    FFT.C
    FourierMap.C
    ScalarFunction.C
    VectorFunction.C
)

set_property(TARGET qcMath PROPERTY CXX_STANDARD 20)
target_compile_options(qcMath PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcMath PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})

target_include_directories(qcMath PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_link_libraries(qcMath PUBLIC)
