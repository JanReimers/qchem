add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
BasisSet.C
TIBS_Common.C
HF_IBS_Common.C
DHF_IBS_Common.C
IBS_Common.C
DFT_IBS_Common.C
Fit_IBS_Common.C
BS_Common.C
Cache4.C
DataBase/HeapDB.C
ERI4.C
Molecule/Factory.C
Molecule/Unit.C
Molecule/PolarizedGaussian/Triangle3D.C
Molecule/PolarizedGaussian/AuxillaryData.C
Molecule/PolarizedGaussian/AuxillaryFJ.C
Molecule/PolarizedGaussian/CDCache.C
Molecule/PolarizedGaussian/IrrepBasisSet.C
Molecule/PolarizedGaussian/Block.C
Molecule/PolarizedGaussian/MnD/Hermite2.C
Molecule/PolarizedGaussian/MnD/Hermite1.C
Molecule/PolarizedGaussian/MnD/RNLM.C
Molecule/PolarizedGaussian/Readers/Gaussian94.C
Molecule/PolarizedGaussian/IEClient.C
Molecule/PolarizedGaussian/IntegralEngine.C
Molecule/PolarizedGaussian/Radial/GaussianH3.C
Molecule/PolarizedGaussian/Radial/RadialFunctionImplementation.C
Molecule/PolarizedGaussian/Radial/ContractedGaussianH3.C
Molecule/PolarizedGaussian/Radial/GaussianCD.C
Molecule/PolarizedGaussian/Radial/ContractedGaussianRF.C
Molecule/PolarizedGaussian/Radial/GaussianRF.C
Molecule/PolarizedGaussian/BasisFunction.C
Molecule/PolarizedGaussian/Polarization.C
Molecule/PolarizedGaussian/BasisSet.C
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

