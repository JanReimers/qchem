// File: PolarizedGaussian/IrrepBasisSet.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussian_IrrepBasisSet_H_
#define _PolarizedGaussian_IrrepBasisSet_H_

#include "Polarization.H"
#include "IntegralEngine.H"
#include "IEClient.H"
#include "Block.H"
#include <vector>
#include <memory>
import qchem.BasisSet.Internal.IBS_Common;
import qchem.DFT_IBS;
import qchem.HF_IBS;

namespace PolarizedGaussian
{
    class Reader;

    // Common implementation for orbital and fit basis sets.
    class IrrepBasisSet
        : public virtual ::IrrepBasisSet,
          public IBS_Common,
          public IrrepIEClient
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Cluster *);
        IrrepBasisSet(const RVec &exponents, size_t L, const Cluster *);
        IrrepBasisSet(const RVec &exponents, size_t L);

        virtual std::ostream &Write(std::ostream &) const;

    private:
        
        IrrepBasisSet(const IrrepBasisSet *bs, const bv_t&);
        void MakeBasisFunctions(const RVec &norms);

        bv_t itsBlocks;
    };

    class Orbital_IBS
        : public virtual TOrbital_HF_IBS<double>,
          public virtual TOrbital_DFT_IBS<double>,
          public IrrepBasisSet,
          public Orbital_IBS_Common<double>,
          public Orbital_DFT_IBS_Common<double>,
          public Orbital_HF_IBS_Common<double>,
          public Orbital_IE

    {
        typedef DB_BS_2E<double> db_t;
    public:
        Orbital_IBS(const db_t* db, Reader *, const Cluster *);
        Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L, const Cluster *);
        Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L);

        virtual ::Fit_IBS *CreateCDFitBasisSet(const ::BasisSet*,const Cluster *) const;
        virtual ::Fit_IBS *CreateVxcFitBasisSet(const ::BasisSet*,const Cluster *) const;
        virtual IrrepBasisSet *Clone(const RVec3 &) const;
    };

    class Fit_IBS
        : public virtual ::Fit_IBS
        , public virtual FitIntegrals
        , public IrrepBasisSet
        , public TIBS_Common<double>
        , public Fit_IBS_Common
        , public Fit_IE
    {
    public:
        Fit_IBS(const DB_cache<double>* db , Reader *, const Cluster *);

        virtual ::Fit_IBS *Clone(const RVec3 &) const;
    };

} // namespace PolarizedGaussian

#endif //_PolarizedGaussian_IrrepBasisSet_H_
