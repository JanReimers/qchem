add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
BasisSet.C
     BS_Common.C
    IBS_Common.C
   TIBS_Common.C
 HF_IBS_Common.C
DHF_IBS_Common.C
DFT_IBS_Common.C
Fit_IBS_Common.C
DataBase/HeapDB.C
Cache4.C
ERI4.C
)

target_include_directories(qcBasisSet PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/BasisSet 
    ${PROJECT_SOURCE_DIR}/submodules/oml/include 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)

add_subdirectory(Atom)
add_subdirectory(Molecule)

