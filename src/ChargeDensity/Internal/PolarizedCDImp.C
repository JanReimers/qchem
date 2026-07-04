// File: PolarizedCDImp.C  Implementation for a polarized charge density.
export module qchem.ChargeDensity.Imp.PolarizedCD;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Spin;

export namespace qchem::ChargeDensity
{

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
class Polarized_CDImp
    : public virtual Polarized_CD
{
public:

    Polarized_CDImp(); // No UT coverage
    Polarized_CDImp(rDM_CD* up,rDM_CD* down);
    ~Polarized_CDImp();

          rDM_CD* GetChargeDensity(const Spin&)      ;
    const rDM_CD* GetChargeDensity(const Spin&) const;
    
private:
    rDM_CD* itsSpinUpCD;
    rDM_CD* itsSpinDownCD;
};

} //namespace