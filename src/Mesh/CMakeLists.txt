add_library(qcMesh) 
target_sources(qcMesh PRIVATE Factory.C EulerMaclarenAngularMesh.C GaussLegendreAngularMesh.C LinearMesh.C Mesh.C RadialMesh.C
GaussAngularMesh.C GaussLegendre.C LogRadialMesh.C  MeshIntegrator.C MHLRadialMesh.C  ScalarFunction.C VectorFunction.C)
target_include_directories(qcMesh PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Mesh
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    )
