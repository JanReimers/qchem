add_library(qc_Mesh) 
target_sources(qc_Mesh PRIVATE EulerMaclarenAngularMesh.C GaussLegendreAngularMesh.C LinearMesh.C Mesh.C RadialMesh.C
GaussAngularMesh.C GaussLegendre.C LogRadialMesh.C  MHLRadialMesh.C)
target_include_directories(qc_Mesh PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
