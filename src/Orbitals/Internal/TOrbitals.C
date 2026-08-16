// File: TOrbitals.C  
module;
#include <limits>
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

using qchem::ChargeDensity::rDM_CD;

template <class T> class TOrbitalsImp
    : public virtual Orbitals
    , public virtual TOrbitals<T>
{
    typedef typename TOrbitals<T>::ds_t ds_t; //{double,smat_t}
public:
    TOrbitalsImp(const tobs_t<T>*, Spin s);
    virtual ~TOrbitalsImp();


    virtual ds_t      TakeElectrons      (double ne      )      ;
    virtual ds_t      TakeElectrons      (double ne, const rvec_t& priority);
    virtual ds_t      TakeElectronsFermi (double ne, double kT);
    virtual ds_t      TakeElectronsFermi (double ne, double kT, const rvec_t& eShift);
    virtual ds_t      SetFermiOccupationsAtMu (double mu, double kT, const rvec_t& eShift);
    virtual double    GetChemicalPotential   () const {return itsMu;}
    virtual size_t    GetNumOrbitals     (               ) const;
    virtual size_t    GetNumOccOrbitals  (               ) const;
    virtual double    GetEigenValueChange(const Orbitals&) const;
    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity() const;   //!< BUILDS it (V1.25: owning return)
    virtual rvec_t    GetBasisPopulations(const hmat_t<T>& S) const;
    virtual void      UpdateOrbitals     (const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e);
    virtual Irrep GetQNs() const;
    
    virtual vec_t    <T> operator() (const rvec3_t&) const;
    virtual vec3vec_t<T> Gradient   (const rvec3_t&) const;


    virtual const Orbital* GetOrbital(size_t i) const {return itsOrbitals[i].get();}
    virtual       Orbital* GetOrbital(size_t i)       {return itsOrbitals[i].get();}


    virtual std::ostream&          Write(std::ostream&) const;

private:
    TOrbitalsImp(const TOrbitalsImp&);
    ds_t BuildDensity(double ne);   //!< build D/D' from the occupied orbitals (shared by both TakeElectrons)

    //! μ of the last Fermi fill -- captured in SetFermiOccupationsAtMu, the ONE funnel both the per-block
    //! solve and the shared-reservoir fill pass through.  NaN until a Fermi fill runs (a cold fill has no μ).
    double            itsMu=std::numeric_limits<double>::quiet_NaN();
    const tobs_t<T>*  itsBasisSet;
    std::vector<std::unique_ptr<Orbital>> itsOrbitals;
    Irrep         itsQNs;
    hmat_t<T>         itsD; // density matrix D=C*Cd (outer product); Hermitian (= symmetric for real T)
    size_t            itsM; // orthonormal-basis dim (columns of the ortho transform V): n, or n-k after a
                            // canonical-ortho truncation.  Sizes the orthonormal D' (D stays n x n).
};

} //namespace