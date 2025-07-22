// File: BFGrouper.C  Group Slater or Gaussian basis functions by unique exponents.
module;
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
module qchem.BasisSet.Atom.BFGrouper;
import qchem.BasisSet.Atom.IEClient;

using std::cout;
using std::endl;

void BFGrouper::Append(AtomIrrepIEClient* aiec)
{
    for (auto e:aiec->es)
    {
        size_t index=unique_es.size();
        if (const auto &ie =unique_es.find(e);ie==unique_es.end())
        {
            unique_esv.push_back(e);
            maxls.push_back(aiec->l);
            unique_es[e]=index;
//            for (auto e:unique_esv) cout << e << " ";
//            cout << endl;
        }
        else 
            index=ie->second;

        if (aiec->l>maxls[index]) maxls[index]=aiec->l;
//        cout << "BFGrouper index,l,maxl=" << index << " " << aiec->l << " " << maxls[index] << endl;
        aiec->es_indices.push_back(index);  
    }        
}
 


//  indices should already be zero based.
size_t BFGrouper::LMax(size_t ia, size_t ib, size_t ic, size_t id) const
{
    size_t lmax_ab=std::max(maxls[ia],maxls[ib]);
    size_t lmax_cd=std::max(maxls[ic],maxls[id]);
//    cout << "(abcd)=(" << ia << " "  << ib << " "  << ic << " "  << id << ") max(ab)=" << lmax_ab << " max(cd)=" << lmax_cd <<endl;
    return std::max(lmax_ab,lmax_cd);
}


