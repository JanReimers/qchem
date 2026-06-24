// File: TOrbitals.C  
module;
#include <iosfwd>
#include <vector>
#include <memory>
export module qchem.Orbitals.Internal.OrbitalsImp;
export import qchem.Orbitals;
export import qchem.Symmetry.Irrep;
export import qchem.Types;

import qchem.Orbitals.Types;


export namespace qchem::Orbitals
{

using qchem::ChargeDensity::DM_CD;

template <class T> class TOrbitalsImp
    : public virtual Orbitals
    , public virtual TOrbitals<T>
{
    typedef typename TOrbitals<T>::ds_t ds_t; //{double,smat_t}
public:
    TOrbitalsImp(const tobs_t<T>*, Spin s);
    virtual ~TOrbitalsImp();


    virtual ds_t      TakeElectrons      (double ne      )      ;
    virtual size_t    GetNumOrbitals     (               ) const;
    virtual size_t    GetNumOccOrbitals  (               ) const;
    virtual double    GetEigenValueChange(const Orbitals&) const;
    virtual tDM_CD<T>* GetChargeDensity  () const;
    virtual void      UpdateOrbitals     (const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e);
    virtual Irrep GetQNs() const;
    
    virtual vec_t    <T> operator() (const rvec3_t&) const;
    virtual vec3vec_t<T> Gradient   (const rvec3_t&) const;


    virtual const Orbital* GetOrbital(size_t i) const {return itsOrbitals[i].get();}
    virtual       Orbital* GetOrbital(size_t i)       {return itsOrbitals[i].get();}


    virtual std::ostream&          Write(std::ostream&) const;

private:
    TOrbitalsImp(const TOrbitalsImp&);

    const tobs_t<T>*  itsBasisSet;
    std::vector<std::unique_ptr<Orbital>> itsOrbitals;
    Irrep         itsQNs;
    hmat_t<T>         itsD; // density matrix D=C*Cd (outer product); Hermitian (= symmetric for real T)
};

} //namespace