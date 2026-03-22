add_library(qcLASolver_blaze) 
target_sources(qcLASolver_blaze PRIVATE 
    Internal/Imp/LASolverImp.C  
    Internal/Imp/LASolverLapack.C  
    Imp/Factory.C
    )
target_sources(qcLASolver_blaze PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    LASolver.C
    Internal/LASolverImp.C  
    Internal/LASolverLapack.C  
)
set_property(TARGET qcLASolver_blaze PROPERTY CXX_STANDARD 20)
target_include_directories(qcLASolver_blaze PRIVATE 
      ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_link_libraries(qcLASolver_blaze PUBLIC qcLASolver qcCommon) 

