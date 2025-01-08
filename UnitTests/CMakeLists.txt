cmake_minimum_required(VERSION 3.5.0)

file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
link_directories(../../oml/Debug)
add_executable(UTMain UTMain.C DiracIntegral.C A_DHF.C A_HF.C QchemTester.C)

target_link_options(UTMain PRIVATE -fopenmp)
target_include_directories(UTMain PRIVATE ../include ../../oml/include)
target_link_libraries(UTMain qchem)
