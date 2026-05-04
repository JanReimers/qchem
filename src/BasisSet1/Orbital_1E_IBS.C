// File: BasisSet/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <cassert>
#include <string>
export module qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.IrrepBasisSet;
export import qchem.Cluster;

import qchem.BasisSet1.DB_Cache;

export namespace BasisSet1
{
//  The are used for caching 1) radial Slater integrals R_k(abcd) 2) Direct/Exchange integrals
class IrrepBasisSet_IDs
{
public:
    virtual std::string  RadialID() const=0;
    virtual std::string AngularID() const=0;
    virtual std::string Name     () const=0;

};

//--------------------------------------------------------------------------------
//
//! The method of integral evaluation is of course strongly dependant on the
//! precise details of basis functions or basis set.  
//! All integral functions except MakeNormalization return normalized integrals.
//! Interfaces for 1 electron integrals used for all Irrep basis sets: 1E,Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for overlap integrals.
//! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
template <class T> class Integrals_Overlap : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T>  MakeOverlap() const=0; //Only called once for a given {radial,angular} ID pair.
    const   smat_t<T>&     Overlap() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
            ? cache->GetSMat() : cache->Set(MakeOverlap());
    }
};

//! \brief Interface for Kinetic energy integrals.
//! Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
template <class T> class Integrals_Kinetic : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeKinetic() const=0; //Only called once for a given {radial,angular} ID pair.
    
    const smat_t<T>& Kinetic() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2C::Kinetic,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
            ? cache->GetSMat() : cache->Set(MakeKinetic());
       
    }
};

//! \brief Interface for electron-nucleus attraction integrals.
//! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
template <class T> class Integrals_Nuclear : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeNuclear(const Cluster* cl) const=0; //Only called once for a given {radial,angular, cluster} ID triple.
    
    const smat_t<T>& Nuclear(const Cluster* cl) const
    {
        assert(cl);
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2n::Nuclear,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),cl->ID())
            ? cache->GetSMat() : cache->Set(MakeNuclear(cl));
    }
};


//
// Define an orbital irrep basis set which supports integrals for 1-electron (1E) orbital calculations.
// Mix-in the integral interfaces required for a 1E orbital basis. 
//
template <class T> class Orbital_1E_IBS
    : public IrrepBasisSet<T> //brings in symmetry and op()(r)
    , public virtual Integrals_Overlap<T> 
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
public:    
    Orbital_1E_IBS(const Irrep_QNs::sym_t& sym) : IrrepBasisSet<T>(sym) {};
};

typedef Orbital_1E_IBS<double>    Real_OIBS;
typedef Orbital_1E_IBS<dcmplx> Complex_OIBS;

} //namespace
