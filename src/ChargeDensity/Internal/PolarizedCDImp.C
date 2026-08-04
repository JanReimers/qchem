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
//  Templated like tPolarized_CD: <double> = molecular, <dcmplx> = polarized plane-wave (tier 4b).
//
template <class T> class tPolarized_CDImp
    : public virtual tPolarized_CD<T>
{
public:

    tPolarized_CDImp(); // No UT coverage
    tPolarized_CDImp(tDM_CD<T>* up,tDM_CD<T>* down);
    ~tPolarized_CDImp();

          tDM_CD<T>* GetChargeDensity(const Spin&)      ;
    const tDM_CD<T>* GetChargeDensity(const Spin&) const;

private:
    tDM_CD<T>* itsSpinUpCD;
    tDM_CD<T>* itsSpinDownCD;
};

using Polarized_CDImp = tPolarized_CDImp<double>;

} //namespace
