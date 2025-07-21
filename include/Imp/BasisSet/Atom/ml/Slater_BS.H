// File: Atom/ml/Slater_BS.H  r^l exp(-a*r)*Y_lm type basis set.
#ifndef _Atom_ml_Slater_BS_H_
#define _Atom_ml_Slater_BS_H_

#include "../l/Slater_BS.H"
import qchem.BasisSet.Atom.ml.Angular;
import qchem.BasisSet.Common;
import qchem.Symmetry.ElectronConfiguration;

namespace Atom_ml
{
namespace Slater
{

class BasisSet 
    : public ::Slater::BS_Common
    , public IE_BS_2E_Angular
{
public:
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
    
};

}} //namespace Slater_m

#endif // _Atom_ml_Slater_BS_H_
