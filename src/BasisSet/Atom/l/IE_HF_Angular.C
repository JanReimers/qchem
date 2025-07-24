// File: Atoml_IE_HF_Angular.C  Angular 2e-Integrals for atoml HF basis sets.
module;
#include <cmath>
export module qchem.BasisSet.Atom.Internal.l.Angular;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.Internal.AngularIntegrals;

export namespace Atoml
{

class IE_BS_2E_Angular : public virtual ::AtomIE_BS_2E_Angular
{
public:
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const
    {
        return AngularIntegrals::Coulomb(a->l,c->l);
    }
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const
    {
        return AngularIntegrals::Exchange(a->l,b->l);
    }

};

 

}