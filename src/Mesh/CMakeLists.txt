add_library(qcMesh) 
target_sources(qcMesh PRIVATE 
    Imp/Factory.C   
    Imp/ScalarFunction.C 
    Imp/VectorFunction.C 
    Imp/MeshIntegrator.C 
    Internal/LinearMesh.C 
    Internal/GaussAngularMesh.C 
    Internal/GaussLegendreAngularMesh.C 
    Internal/EulerMaclarenAngularMesh.C 
    Internal/MHLRadialMesh.C 
    Internal/LogRadialMesh.C 
    )
target_sources(qcMesh PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Mesh.C 
    RadialMesh.C 
    MeshIntegrator.C 
    ScalarFunction.C
    VectorFunction.C
    Factory.C
    Internal/MeshTypes.C 
    Internal/RadialMeshTypes.C  
     )
target_include_directories(qcMesh PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    )
# target_compile_options(qcMesh PRIVATE -Wall -Wno-unknown-pragmas -fmodules-ts
#     $<$<CONFIG:Debug>:-g -O0>
#     $<$<CONFIG:Release>: -O2>
# )
# We have to be on cxx-20 or high again for cmake/gcc to find .gcm binaru module  interfact (BMI) files!
# If you leave this out it does not complain about the import statement, in compains that it can;t find the .gcm file.
set_property(TARGET qcMesh PROPERTY CXX_STANDARD 20)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
target_link_libraries(qcMesh PUBLIC qcCommon omlm) 
