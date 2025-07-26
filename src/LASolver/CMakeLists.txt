add_library(qcLASolver) 
target_sources(qcLASolver PRIVATE Imp/LAParams.C  Internal/Imp/LASolverImp.C  Internal/Imp/LASolverLapack.C  Internal/Imp/LASolverOML.C Imp/Factory.C)
target_sources(qcLASolver PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    LAParams.C
    LASolver.C
    Internal/LASolverImp.C  
    Internal/LASolverLapack.C  
    Internal/LASolverOML.C
)
set_property(TARGET qcLASolver PROPERTY CXX_STANDARD 20)
target_link_libraries(qcLASolver PUBLIC qcCommon omlm) 

