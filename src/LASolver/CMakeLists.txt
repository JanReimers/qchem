add_library(qcLASolver) 
target_sources(qcLASolver PRIVATE 
    Internal/Imp/LASolverImp.C  
    Internal/Imp/LASolverLapack.C  
    Imp/Factory.C
    )
target_sources(qcLASolver PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    LASolver.C
    LAParams.C
    Internal/LASolverImp.C  
    Internal/LASolverLapack.C  
)
set_property(TARGET qcLASolver PROPERTY CXX_STANDARD 20)
target_include_directories(qcLASolver PRIVATE 
      ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_link_libraries(qcLASolver PUBLIC qcCommon) 

