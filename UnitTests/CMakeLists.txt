cmake_minimum_required(VERSION 3.5.0)

file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../../oml/Debug>
    $<$<CONFIG:Release>: ../../oml/Release>
)
add_executable(UTMain UTMain.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C M_PG_U.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C QchemTester.C)
#add_executable(UTMain UTMain.C BSplines.C QchemTester.C)

target_link_options(UTMain PRIVATE -fopenmp -pg)
target_include_directories(UTMain PRIVATE ../include ../../oml/include ../BSplinebasis/include)
target_link_libraries(UTMain qchem)
