// File:: Symmetry/Imp/Factory.C  Create instances for various symmetry classes.
module;
#include <vector>
module qchem.Symmetry.Factory;
// We want to try hide all these imports from the main code.
import qchem.Symmetry.Internal.Spherical;
import qchem.Symmetry.BlochQN;
import qchem.Symmetry.Unit;

using namespace Symmetryns::Internal::Spherical;

namespace SymmetryFactory
{
sym_t YFactory(size_t l,const ivec_t& mls)
{
    if (mls.size()==0)
        return sym_t(new Yl_Sym(l));
    else
        return sym_t(new Ylm_Sym(l,mls));
}

sym_t     ΩFactory(int κ,const rvec_t& mjs)
{
     if (mjs.size()==0)
        return sym_t(new Omega_k_Sym(κ));
    else
        return sym_t(new Omega_kmj_Sym(κ,mjs));
}
sym_t BlochFactory(ivec3_t N, ivec3_t k)
{
    return sym_t(new BlochQN(N,k));
}
sym_t  UnitFactory()
{
    return sym_t(new UnitQN);
}
    
} //namesapce
