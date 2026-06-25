add_library(qcMesh)
target_sources(qcMesh PRIVATE
    Imp/ScalarFunction.C
    Imp/VectorFunction.C
    )
target_sources(qcMesh PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Mesh.C
    ScalarFunction.C
    VectorFunction.C
     )
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
target_include_directories(qcMesh PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# This linker association is required in order for cmake/gcc to find .gcm binary module interface (BMI) files.
target_link_libraries(qcMesh PUBLIC qcCommon)
