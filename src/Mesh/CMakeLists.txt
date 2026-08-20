add_library(qcMesh)
# Module implementation units (no module interface; not part of the BMI file set).
target_sources(qcMesh PRIVATE
    Internal/MHLRadialMesh.C
    Internal/LogRadialMesh.C
    Internal/LinearRadialMesh.C
    Internal/RadialFactory.C
    Internal/LebedevAngularMesh.C
    Internal/GaussLegendreAngularMesh.C
    Internal/AngularFactory.C
    )
# Module interface units (export the BMI .gcm files other modules import).
target_sources(qcMesh PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Mesh.C
    GaussLegendre.C
    RadialMesh.C
    AngularMesh.C
    Product.C
    Quadrature.C
    Builder.C
    XCPolicy.C
    )
# Flag uniformity is a BUILD CORRECTNESS requirement under C++20 modules: a .pcm and every TU that
# consumes it must agree on target features, so a lib that opts out of the common set makes any
# global ISA flag unbuildable (measured 2026-08-20).  This lib previously set NO options at all and
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
target_include_directories(qcMesh PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# Linking qcCommon is also what lets cmake/gcc find its exported .gcm (BMI) files.
target_link_libraries(qcMesh PUBLIC qcMath qcCommon)

add_subdirectory(tests)
