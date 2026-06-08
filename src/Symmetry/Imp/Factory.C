// File:: Symmetry/Imp/Factory.C  Create instances for various symmetry classes.
module;
#include <vector>
module qchem.Symmetry.Factory;
// We want to try hide all these imports from the main code.
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;
import qchem.Symmetry.Okmj;
import qchem.Symmetry.BlochQN;
import qchem.Symmetry.Unit;
//  import qchem.Symmetry;

namespace SymmetryFactory
{
sym_t YFactory(size_t l,std::vector<int> mls)
{
    if (mls.size()==0)
        return sym_t(new Yl_Sym(l));
    else
        return sym_t(new Ylm_Sym(l,mls));
}

sym_t     ΩFactory(int κ,std::vector<double> mjs)
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
