add_library(qcMesh)
# Module implementation units (no module interface; not part of the BMI file set).
target_sources(qcMesh PRIVATE
    Internal/MHLRadialMesh.C
    Internal/LogRadialMesh.C
    Internal/LinearRadialMesh.C
    Internal/RadialFactory.C
    Internal/GaussAngularMesh.C
    Internal/GaussLegendreAngularMesh.C
    Internal/EulerMaclarenAngularMesh.C
    Internal/AngularFactory.C
    )
# Module interface units (export the BMI .gcm files other modules import).
target_sources(qcMesh PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Mesh.C
    Fields.C
    GaussLegendre.C
    RadialMesh.C
    AngularMesh.C
    Product.C
    Quadrature.C
    Builder.C
    )
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
target_include_directories(qcMesh PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# Linking qcCommon is also what lets cmake/gcc find its exported .gcm (BMI) files.
target_link_libraries(qcMesh PUBLIC qcCommon)

add_subdirectory(tests)
