// File: BSpline/BFGrouper.C  Group BSpline basis functions by support start positions.
module;
#include <cassert>
#include <iostream>
module qchem.Basisset.Atom.radial.BSpline.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;
import qchem.stl_io;
import qchem.Basisset.Atom.radial.BSpline.IEC;


using std::cout;
using std::endl;

namespace BSpline
{
template <size_t K> void BFGrouper<K>::Append(IrrepIEClient<K>* iec)
{
    for (auto s:iec->splines)
    {
        size_t index=unique_sp.size();
        if (const auto &ie =unique_sp.find(s);ie==unique_sp.end())
        {
            unique_spv.push_back(s);
            maxls.push_back(iec->l);
            unique_sp[s]=index;
//            for (auto e:unique_esv) cout << e << " ";
//            cout << endl;
        }
        else 
            index=ie->second;

        if (iec->l>maxls[index]) maxls[index]=iec->l;
//        cout << "BFGrouper index,l,maxl=" << index << " " << aiec->l << " " << maxls[index] << endl;
        iec->sp_indices.push_back(index);  
    }     
    // cout << "l, iec->sp_indices = " << iec->l << " " << iec->sp_indices << endl;
    itsGLs[iec->l]=iec->itsGL;
}
//  indices should already be zero based.
template <size_t K> size_t BFGrouper<K>::LMax(size_t ia, size_t ib, size_t ic, size_t id) const
{
    size_t lmax_ab=std::max(maxls[ia],maxls[ib]);
    size_t lmax_cd=std::max(maxls[ic],maxls[id]);
//    cout << "(abcd)=(" << ia << " "  << ib << " "  << ic << " "  << id << ") max(ab)=" << lmax_ab << " max(cd)=" << lmax_cd <<endl;
    return std::max(lmax_ab,lmax_cd);
}
template <size_t K> const GLCache* BFGrouper<K>::GetGL(size_t l) const
{
    auto i=itsGLs.find(l);
    assert(i!=itsGLs.end());
    return i->second;
}

#define INSTANCEk(k) template class BFGrouper<k>;
#include "../Instance.hpp"

} // namespace
