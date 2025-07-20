// File: PolarizedCDImp.C  Implementation for a polarized charge density.
export module qchem.ChargeDensity.Imp.PolarizedCD;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Spin;
//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
export class Polarized_CDImp
    : public virtual Polarized_CD
{
public:

    Polarized_CDImp(); // No UT coverage
    Polarized_CDImp(DM_CD* up,DM_CD* down);
    ~Polarized_CDImp();

          DM_CD* GetChargeDensity(const Spin&)      ;
    const DM_CD* GetChargeDensity(const Spin&) const;
    
private:
    DM_CD* itsSpinUpCD;
    DM_CD* itsSpinDownCD;
};

