// File: Atom/ml/BSpline_BS.H BSpline Basis Set for atoms, no m degeneracy.
#ifndef _Atom_ml_BSpline_BS_H_
#define _Atom_ml_BSpline_BS_H_

#include "IE_HF_Angular.H"
#include "../../BS_Common.H"
#include "../l/Slater_BS.H"

namespace Atom_ml
{
namespace BSpline
{

template <size_t K> class BasisSet 
    : public ::BSpline::BS_Common<K>
    , public IE_BS_2E_Angular
{
public:
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
    
};

}} //namespace Atom_ml::BSpline

#endif // _Atom_ml_BSpline_BS_H_
