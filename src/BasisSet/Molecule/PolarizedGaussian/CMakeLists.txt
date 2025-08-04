target_sources(qcMolecule_BS PRIVATE 
Imp/BasisSet.C
Imp/IrrepBasisSet.C
MnD/Imp/Triangle3D.C
MnD/Imp/AuxillaryData.C
MnD/Imp/AuxillaryFJ.C
Internal/Imp/CDCache.C
Internal/Imp/Block.C
Internal/Imp/IEClient.C
Internal/Imp/IntegralEngine.C
Internal/Imp/Polarization.C
MnD/Imp/Hermite2.C
MnD/Imp/Hermite1.C
MnD/Imp/RNLM.C
Readers/Imp/Gaussian94.C
Radial/Imp/GaussianH3.C
Radial/Imp/Common.C
Radial/Imp/ContractedGaussianH3.C
Radial/GaussianCD.C
Radial/Imp/ContractedGaussianRF.C
Radial/Imp/GaussianRF.C


)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    BasisSet.C
    RadialFunction.C;
    Reader.C

    MnD/AuxillaryFJ.C
    Internal/CDCache.C
    Internal/Block.C
    Internal/IEClient.C
    Internal/IntegralEngine.C
    Internal/Polarization.C
    MnD/Triangle.C
    MnD/Hermite1.C
    MnD/Hermite2.C
    MnD/Hermite3.C
    MnD/RNLM.C

    Radial/GData.C
    Radial/Common.C
    Radial/GaussianH3.C
    Radial/ContractedGaussianH3.C
    Radial/GaussianRF.C
    Radial/ContractedGaussianRF.C
    Readers/Gaussian94.C
)
 