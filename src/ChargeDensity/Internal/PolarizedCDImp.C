// File: PolarizedCDImp.C  Implementation for a polarized charge density.
module;
#include <memory>   // unique_ptr: this class OWNS its two spin channels (V1.25)
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
    //! TAKES OWNERSHIP of both channels (V1.25).  unique_ptr rather than raw pointers plus a hand-written
    //! dtor: the old form OWNED them (the dtor deleted both) but declared no copy ctor or assignment, so an
    //! implicit copy was a double-delete.  unique_ptr states the ownership and deletes the copy for us.
    tPolarized_CDImp(std::unique_ptr<tDM_CD<T>> up, std::unique_ptr<tDM_CD<T>> down);

          tDM_CD<T>* GetChargeDensity(const Spin&)      ;
    const tDM_CD<T>* GetChargeDensity(const Spin&) const;

private:
    std::unique_ptr<tDM_CD<T>> itsSpinUpCD;
    std::unique_ptr<tDM_CD<T>> itsSpinDownCD;
};

using Polarized_CDImp = tPolarized_CDImp<double>;

} //namespace
