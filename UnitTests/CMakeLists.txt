# cmake_minimum_required(VERSION 3.30.0)

#variable_watch(__CMAKE_CXX_COMPILER_OUTPUT)

#file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>: ../submodules/build/oml/Debug submodules/wignerSymbols/build>
    $<$<CONFIG:Release>: ../submodules/build/oml/Release submodules/wignerSymbols/build>
)

#  Main
add_executable(UTMain UTMain.C json.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C M_PG_U.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C QchemTester.C)
#add_executable(UTMain UTMain.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C QchemTester.C)
#GCM_PATH := $(LIB_PATH)/gcm.cache
#MAPPER := '-fmodule-mapper=|@g++-mapper-server -r'$(GCM_PATH)
#-Wno-unknown-pragmas
# -Wno-cpp  Suppress warning messages emitted by #warning directives.

target_compile_options(UTMain PRIVATE -Wall -Wno-unknown-pragmas -Wno-cpp 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
set_property(TARGET UTMain PROPERTY CXX_STANDARD 20)

target_link_options(UTMain PRIVATE )
target_include_directories(UTMain PRIVATE 
    ../submodules/BSplinebasis/include
      ${PROJECT_SOURCE_DIR}/submodules/json/single_include

)
target_link_libraries(UTMain PUBLIC qchem qcCommon omlm)

# Symmetry module
add_executable(UTSym UTSym.C SymQN.C ElectronConfigurations.C)
target_link_libraries(UTSym PUBLIC qcSymmetry gtest)
set_property(TARGET UTSym PROPERTY CXX_STANDARD 20)
# Cluster module
add_executable(UTCluster UTCluster.C)
target_link_libraries(UTCluster PUBLIC qcCluster gtest)
set_property(TARGET UTCluster PROPERTY CXX_STANDARD 20)
# BasisSet module
add_executable(UTBasisSet UTBasisSet.C BasisSet_SL.C)
target_link_libraries(UTBasisSet PUBLIC qcAtom_BS wignerSymbols gtest) #qcLattice_BS
set_property(TARGET UTBasisSet PROPERTY CXX_STANDARD 20)

target_compile_definitions(UTBasisSet PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_compile_options(UTBasisSet PRIVATE -Wall -Wno-unknown-pragmas -Wno-cpp 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)



# setup testing
enable_testing()
include(CTest)
include(GoogleTest)
add_test(NAME Symmetry COMMAND UTSym )
gtest_discover_tests(UTSym) 