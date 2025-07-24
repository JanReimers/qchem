target_sources(qcMolecule_BS PRIVATE 
Imp/BasisSet.C
Imp/IrrepBasisSet.C
Imp/BasisFunction.C
Triangle3D.C
MnD/Imp/AuxillaryData.C
MnD/Imp/AuxillaryFJ.C
Internal/Imp/CDCache.C
Internal/Imp/Block.C
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
    MnD/AuxillaryFJ.C
    Internal/CDCache.C
    Internal/Block.C
)
 