add_library(qcMesh) 
target_sources(qcMesh PRIVATE Factory.C EulerMaclarenAngularMesh.C GaussLegendreAngularMesh.C LinearMesh.C Mesh.C RadialMesh.C
GaussAngularMesh.C GaussLegendre.C LogRadialMesh.C  MeshIntegrator.C MHLRadialMesh.C  ScalarFunction.C VectorFunction.C)
target_include_directories(qcMesh PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Mesh
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    )
target_compile_options(qcMesh PRIVATE -Wall -Wno-unknown-pragmas -fmodules
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
