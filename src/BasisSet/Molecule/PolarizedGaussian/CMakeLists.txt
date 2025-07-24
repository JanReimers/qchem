target_sources(qcMolecule_BS PRIVATE 
Imp/BasisSet.C
Imp/IrrepBasisSet.C
Imp/BasisFunction.C
Triangle3D.C
AuxillaryData.C
AuxillaryFJ.C
CDCache.C
Block.C
MnD/Hermite2.C
MnD/Hermite1.C
MnD/RNLM.C
Readers/Gaussian94.C
IEClient.C
IntegralEngine.C
Radial/GaussianH3.C
Radial/RadialFunctionImplementation.C
Radial/ContractedGaussianH3.C
Radial/GaussianCD.C
Radial/ContractedGaussianRF.C
Radial/GaussianRF.C

Polarization.C

)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
   BasisSet.C
)
 