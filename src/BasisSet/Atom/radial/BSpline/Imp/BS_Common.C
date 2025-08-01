// File: Atom/radial/BSpline/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.
module;
#include <cassert>
module qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.BasisSet.Internal.IEClient;
import qchem.Irrep_BS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;
import oml;

namespace BSpline
{
    template <size_t K> BS_Common<K>::~BS_Common()
    {
        delete itsRkCache;
    }
    template <size_t K> void BS_Common<K>::Insert(bs_t* bs)
    {
        ::BS_Common::Insert(bs);
        auto iec=dynamic_cast<const ::IrrepIEClient*>(bs);
        assert(iec);
        this->Append(iec);

    }
    template <size_t K> void BS_Common<K>::BuildCache(size_t lmax)
    {
        itsRkCache=new RkCache<K>(unique_spv,*this->GetGL(lmax),lmax);
    }
    template <size_t K> const Cacheable* BS_Common<K>::Create(size_t ia,size_t ic,size_t ib,size_t id) const
    {
        assert(itsRkCache);
        // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
        size_t lmax=LMax(ia,ib,ic,id);
        const GLCache* gl=this->GetGL(lmax);
        return new BSpline::RkEngine(unique_spv,ia,ib,ic,id,lmax,*gl,*itsRkCache);

    }
    template <size_t K> Vector<double> BS_Common<K>::loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->Coulomb_Rk(la,lc);
    }
    template <size_t K> Vector<double> BS_Common<K>::loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc);
    }
    
#define INSTANCEk(k) template class BS_Common<k>;
#include "../Instance.hpp"
    
}