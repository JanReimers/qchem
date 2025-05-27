cmake_minimum_required(VERSION 3.5.0)

variable_watch(__CMAKE_CXX_COMPILER_OUTPUT)

file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../../oml/Debug>
    $<$<CONFIG:Release>: ../../oml/Release>
)
add_executable(UTMain UTMain.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C M_PG_U.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C QchemTester.C)
#add_executable(UTMain UTMain.C OpenMP.C)

target_compile_options(UTMain PRIVATE -fopenmp -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UTMain PRIVATE -fopenmp)
target_include_directories(UTMain PRIVATE ../include ../../oml/include ../BSplinebasis/include)
target_link_libraries(UTMain qchem)
