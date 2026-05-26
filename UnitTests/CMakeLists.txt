# cmake_minimum_required(VERSION 3.30.0)

#variable_watch(__CMAKE_CXX_COMPILER_OUTPUT)

#file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories( 
    $<$<CONFIG:Debug>:  submodules/wignerSymbols/build>
    $<$<CONFIG:Release>:  submodules/wignerSymbols/build>
)

#  Main
# 
add_executable(UTMain gtestmain.C json.C BSplines.C EigenSolverUT.C SymQN.C ElectronConfigurations.C 
DiracIntegral.C SlaterIntegral.C GaussianIntegral.C A_DHF.C A_HF.C A_DFT.C M_PG_U.C QchemTester.C Angular.C ERI4.C MeshIntegrals.C
UTCluster.C UTBasisSet.C BasisSet1_Atom.C BasisSet_Atom.C PeriodicTable.C)
#add_executable(UTMain UTMain.C DiracIntegral.C SlaterIntegral.C GaussianIntegral.C QchemTester.C)
#GCM_PATH := $(LIB_PATH)/gcm.cache
#MAPPER := '-fmodule-mapper=|@g++-mapper-server -r'$(GCM_PATH)
#-Wno-unknown-pragmas
# -Wno-cpp  Suppress warning messages emitted by #warning directives.

target_compile_options(UTMain PRIVATE -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
target_compile_definitions(UTMain PRIVATE ${CXX_COMMON_DEFINES}
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
set_property(TARGET UTMain PROPERTY CXX_STANDARD 20)

target_link_options(UTMain PRIVATE -pg)
target_include_directories(UTMain PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include

)
target_link_libraries(UTMain PUBLIC qchem qcCommon)

# -isystem /usr/include/c++/v1
# --gcc-toolchain=/opt/gcc-15
add_executable(UTMain1 gtestmain.C QchemTesterImp.C A_Pool.C BasisSetPoolImp.C)
target_sources(UTMain1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES QchemTester.C BasisSetPool.C)

target_compile_options(UTMain1 PRIVATE -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
target_compile_definitions(UTMain1 PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
set_property(TARGET UTMain1 PROPERTY CXX_STANDARD 20)

target_link_options(UTMain1 PRIVATE -pg )
target_include_directories(UTMain1 PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include

)
target_link_libraries(UTMain1 PUBLIC qchem qcCommon)


# Symmetry module
add_executable(UTSym gtestmain.C SymQN.C ElectronConfigurations.C)
target_include_directories(UTSym PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include

)
target_link_libraries(UTSym PUBLIC qcFactory qcBasisSet1 qcSymmetry gtest lapack wignerSymbols)
set_property(TARGET UTSym PROPERTY CXX_STANDARD 20)
# Cluster module
add_executable(UTCluster gtestmain.C UTCluster.C)
target_link_libraries(UTCluster PUBLIC qcCluster gtest)
set_property(TARGET UTCluster PROPERTY CXX_STANDARD 20)
# BasisSet module
add_executable(UTBasisSet gtestmain.C UTBasisSet.C BasisSet_Atom.C)
target_link_libraries(UTBasisSet PUBLIC qcAtom_BS1 lapack wignerSymbols gtest) #qcLattice_BS
target_include_directories(UTBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
set_property(TARGET UTBasisSet PROPERTY CXX_STANDARD 20)

target_compile_definitions(UTBasisSet PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_compile_options(UTBasisSet PRIVATE -Wall -Wno-unknown-pragmas -Wno-cpp 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

add_executable(HydrogenAI2 HydrogenAI2.C)
target_include_directories(HydrogenAI2 PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(HydrogenAI2 PUBLIC lapack blas)

# setup testing
enable_testing()
include(CTest)
include(GoogleTest)
add_test(NAME Symmetry COMMAND UTSym )
gtest_discover_tests(UTSym) 