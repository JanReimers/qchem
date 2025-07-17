add_library(qcLASolver) 
target_sources(qcLASolver PRIVATE LAParams.C  LASolverImp.C  LASolverLapack.C  LASolverOML.C)
target_include_directories(qcLASolver PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/LASolver  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include )

set_property(TARGET qcLASolver PROPERTY CXX_STANDARD 20)
target_link_libraries(qcLASolver PUBLIC omlm) 

