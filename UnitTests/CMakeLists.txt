cmake_minimum_required(VERSION 3.5.0)

variable_watch(__CMAKE_CXX_COMPILER_OUTPUT)

file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../submodules/build/oml/Debug submodules/wignerSymbols/build>
    $<$<CONFIG:Release>: ../submodules/build/oml/Release submodules/wignerSymbols/build>
)
add_executable(UTMain UTMain.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C M_PG_U.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C QchemTester.C)
#add_executable(UTMain UTMain.C Angular.C)

#-Wno-unknown-pragmas
target_compile_options(UTMain PRIVATE -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UTMain PRIVATE -fopenmp)
target_include_directories(UTMain PRIVATE ../include ../submodules/oml/include ../submodules/BSplinebasis/include)
target_link_libraries(UTMain qchem qc_Common)
