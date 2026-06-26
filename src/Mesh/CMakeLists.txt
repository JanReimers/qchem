add_library(qcMesh)
target_sources(qcMesh PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Mesh.C
     )
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
target_include_directories(qcMesh PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcMesh PUBLIC qcCommon)
