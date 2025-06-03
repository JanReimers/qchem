// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"

SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const SMat& _S) : S(_S) {};
SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() {};
void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las)
{
    itsLaSolver=las;
}

#include <Irrep_BS.H>

SCFAccelerator_DIIS::SCFAccelerator_DIIS() {};
SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const TOrbital_IBS<double>* bs) const
{
    return new SCFIrrepAccelerator_DIIS(bs->Overlap());
}

