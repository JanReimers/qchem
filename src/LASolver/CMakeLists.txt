add_library(qcLASolver) 
target_sources(qcLASolver PRIVATE LAParams.C  LASolverImp.C  LASolverLapack.C  LASolverOML.C)
target_include_directories(qcLASolver PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include )
