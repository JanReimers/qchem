cmake_minimum_required(VERSION 3.5.0)

variable_watch(__CMAKE_CXX_COMPILER_OUTPUT)

#file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../submodules/build/oml/Debug submodules/wignerSymbols/build>
    $<$<CONFIG:Release>: ../submodules/build/oml/Release submodules/wignerSymbols/build>
)

#  Main
add_executable(UTMain UTMain.C json.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C M_PG_U.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C QchemTester.C)

#-Wno-unknown-pragmas
target_compile_options(UTMain PRIVATE -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UTMain PRIVATE -fopenmp)
target_include_directories(UTMain PRIVATE 
    ../include 
    ../submodules/oml/include 
    ../submodules/BSplinebasis/include
      ${PROJECT_SOURCE_DIR}/submodules/json/single_include

)
target_link_libraries(UTMain qchem )

add_executable(UTModules UTModules.C)
set_property(TARGET UTModules PROPERTY CXX_STANDARD 20)
target_compile_options(UTModules PRIVATE -Wall -Wno-unknown-pragmas -fmodules 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
# Symmetry module
add_executable(UTSym UTSym.C SymQN.C ElectronConfigurations.C)
target_include_directories(UTSym PRIVATE 
    ../include 
    ../submodules/oml/include 
)
target_link_libraries(UTSym qcSymmetry gtest)

# setup testing
enable_testing()
include(CTest)
include(GoogleTest)
add_test(NAME Symmetry COMMAND UTSym )
gtest_discover_tests(UTSym) 