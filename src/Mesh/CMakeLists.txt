add_library(qcMesh) 
target_sources(qcMesh PRIVATE EulerMaclarenAngularMesh.C GaussLegendreAngularMesh.C LinearMesh.C Mesh.C RadialMesh.C
GaussAngularMesh.C GaussLegendre.C LogRadialMesh.C  MHLRadialMesh.C)
target_include_directories(qcMesh PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
