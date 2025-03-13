// File: BasisSet/Dirac_IBS.H  Interface for Dirac basis sets with Restricted Kinetic Balance (RKB).

#include "Imp/BasisSet/Dirac_IBS.H"
#include "Imp/Misc/DFTDefines.H"
#include "Imp/Symmetry/OkmjQN.H"

namespace Dirac
{
    template <class T> IntegralEngine<T>
    ::IntegralEngine(const ::Orbital_RKBL_IBS<T>* rkbl, const ::Orbital_RKBS_IBS<T>* rkbs)
    : itsRKBL(rkbl), itsRKBS(rkbs)
    {
        assert(itsRKBL);
        assert(itsRKBS);
    }

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::merge_diag(const SMat& l,const SMat& s)
{
    size_t Nl=l.GetNumRows();
    size_t Ns=s.GetNumRows();
    SMat ls(Nl+Ns);
    Fill(ls,0.0);
    for (auto i:l.rows())
        for (auto j:l.cols(i))
            ls(i,j)=l(i,j);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            ls(Nl+i,Nl+j)=s(i,j);
    return ls;
}

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::merge_off_diag(const Mat& ls)

{
    size_t Nl=ls.GetNumRows();
    size_t Ns=ls.GetNumCols();
    assert(Nl==Ns);
    SMat k(Nl+Ns);
    Fill(k,0.0);
    for (auto i:ls.rows())
        for (auto j:ls.cols())
            k(i,Ns+j)=ls(i,j);
   
    return k;
}    

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::MakeOverlap() const
{
    SMat ol=itsRKBL->Overlap();
    SMat os=itsRKBS->Overlap();
    return merge_diag(ol,os);
}

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::MakeKinetic() const
{
    Mat kls=-2.0*itsRKBL->Kinetic(itsRKBS);
    return merge_off_diag(kls);
}

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::MakeNuclear(const Cluster* c) const
{
    SMat nl=itsRKBL->Nuclear(c);
    SMat ns=itsRKBS->Nuclear(c);
    return merge_diag(nl,ns);
}

template <class T> typename IntegralEngine<T>::SMat IntegralEngine<T>::MakeRestMass() const
{
    static const double f=-2.0*c_light*c_light;
    SMat rl(itsRKBL->size());
    Fill(rl,0.0);
    SMat rs=f*itsRKBS->Overlap();
    return merge_diag(rl,rs);
}


template <class T> ERI4 IntegralEngine<T>::MakeDirect  (const obs_t& c) const
{
    return ERI4();
    // auto da=this;
    // auto dc=dynamic_cast<const Dirac::IntegralEngine<T>*>(&c);
    // assert(da->itsRKBL->GetKappa()==da->itsRKBS->GetKappa());
    // assert(dc->itsRKBL->GetKappa()==dc->itsRKBS->GetKappa());
    // ERI4 JLLLL=itsRKBL->Direct(dc->itsRKBL);
    // ERI4 JLLSS=itsRKBS->Direct(dc->itsRKBS);
    // ERI4 JSSLL=itsRKBS->Direct(dc->itsRKBL);
    // ERI4 JSSSS=itsRKBS->Direct(dc->itsRKBS);
    
    // return merge_diag(JLLLL,JLLSS,JSSLL,JSSSS); 

}
template <class T> ERI4 IntegralEngine<T>::MakeExchange(const obs_t& b) const
{
    return ERI4();
}

template <class T> IrrepBasisSet<T>::IrrepBasisSet(const LAParams& lap
    ,::Orbital_RKBL_IBS<T>* rkbl, int kappa)
    : IrrepBasisSetCommon(new Omega_kQN(kappa))
    , Orbital_IBS_Common<T>(lap)
    , itsRKBL(rkbl)
    , itsRKBS(0)
{
    assert(itsRKBL);
    IntegralEngine<T>::itsRKBL=rkbl;
}
template <class T> IrrepBasisSet<T>::IrrepBasisSet(const LAParams& lap
    ,::Orbital_RKBL_IBS<T>* rkbl,::Orbital_RKBS_IBS<T>* rkbs, int kappa)
    : IrrepBasisSetCommon(new Omega_kQN(kappa))
    , Orbital_IBS_Common<T>(lap)
    , itsRKBL(rkbl)
    , itsRKBS(rkbs)
{
    assert(itsRKBL);
    assert(itsRKBS);
    IntegralEngine<T>::itsRKBL=rkbl;
    IntegralEngine<T>::itsRKBS=rkbs;
}

} // namespace Dirac

template class Dirac::IntegralEngine<double>;
template class Dirac::IrrepBasisSet<double>;