// File: BasisSet/Radial/BSpline/NR/BSpline_BS_Evaluator.C BSpline Basis Set for atoms.
module;
#include <iostream>
#include <map>
#include <cassert>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <type_traits>
export module qchem.BasisSet.Radial.BasisSet;
export import qchem.BasisSet;
export import qchem.ElectronConfiguration;

import qchem.BasisSet.Internal.Orbital_ERI4_IBS;   // the ERI4 substrate the concrete atom IBS realizes
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Radial.Evaluators;
import qchem.Symmetry.Atom.Spherical;   // Symmetry::Atom::Getl (per-l exponent-list ctor)
import qchem.BasisSet.Radial.IBS;
import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Radial.Factory;   // CoulombFitLadder (the face BasisSet_HF realizes)
import qchem.Blaze;

export 
namespace qchem::BasisSet {
namespace Radial {

using namespace Evaluators;

// Full service NR: 1E, DFT and HF
//! The merged, nested Coulomb-fit exponents of \c CoulombFitLadder's rung \a rung >= 1, from every block's fit exponents \a perBlock (already
//! x 2).  Exponents within 2% of a kept neighbour are dropped (a near-duplicate is a near-singular metric, not a richer fit); rung r inserts
//! 2^{r-1} - 1 geometric points between kept neighbours, so each rung's set holds the rung below's.
inline rvec_t MergedCoulombFitExponents(const std::vector<rvec_t>& perBlock, int rung)
{
    std::vector<double> all;
    for (const rvec_t& es : perBlock) for (size_t i=0;i<es.size();i++) all.push_back(es[i]);
    std::sort(all.begin(), all.end());
    std::vector<double> kept;
    for (double e : all) if (kept.empty() || e>kept.back()*1.02) kept.push_back(e);
    const int m=(1<<(rung-1))-1;   // points inserted between neighbours
    std::vector<double> out;
    for (size_t i=0;i<kept.size();i++)
    {
        out.push_back(kept[i]);
        if (i+1<kept.size()) for (int k=1;k<=m;k++) out.push_back(kept[i]*std::pow(kept[i+1]/kept[i], double(k)/(m+1)));
    }
    rvec_t r(out.size());
    for (size_t i=0;i<out.size();i++) r[i]=out[i];
    return r;
}

//! The power-m Coulomb-fit functions of \c CoulombFitLadder's rung \a rung >= 1 for a family whose fit functions carry their own radial
//! power (Slater): each block's fit exponents \a perBlock[b] (already x 2) at half-power \a ls[b] = that block's l, the density's own diagonal
//! products r^{2l} e^{-2 zeta r}; rung r inserts 2^{r-1} - 1 geometric points between each block's neighbours (near-duplicates within a
//! block, 2%, dropped first), so each rung holds the rung below.  Returns (exponents, powers).
inline std::pair<rvec_t,std::vector<int>> ShellPowerCoulombFit(const std::vector<rvec_t>& perBlock, const std::vector<int>& ls, int rung)
{
    std::vector<double> es; std::vector<int> ms;
    const int m=(1<<(rung-1))-1;
    for (size_t b=0;b<perBlock.size();b++)
    {
        std::vector<double> all;
        for (size_t i=0;i<perBlock[b].size();i++) all.push_back(perBlock[b][i]);   // index loop: std iterators cannot see Blaze's op==
        std::sort(all.begin(), all.end());
        std::vector<double> kept;
        for (double e : all) if (kept.empty() || e>kept.back()*1.02) kept.push_back(e);
        for (size_t i=0;i<kept.size();i++)
        {
            es.push_back(kept[i]); ms.push_back(ls[b]);
            if (i+1<kept.size()) for (int k=1;k<=m;k++) { es.push_back(kept[i]*std::pow(kept[i+1]/kept[i], double(k)/(m+1))); ms.push_back(ls[b]); }
        }
    }
    rvec_t r(es.size());
    for (size_t i=0;i<es.size();i++) r[i]=es[i];
    return {r, ms};
}

template <is1E_DFT_HF_Evaluator Evaluator> class BasisSet_HF
    : public virtual Real_BS
    , public BasisSetImp<double>
    , public CoulombFitLadder
{
public:
    // Define the Orbital IBS using mixins.
    class EOrbital_HF_IBS 
        : public Orbital_1E_IBS<Evaluator>
        , public Orbital_DFT_IBS<Evaluator>
        , public Orbital_ERI4_IBS<Evaluator>
        , public IrrepBasisSetImp<Evaluator>
        , public Evaluator
    {
    public:
        EOrbital_HF_IBS(size_t N, double rmin, double rmax, const sym_t& ir)
        : IrrepBasisSetImp<Evaluator>(ir)
        , Evaluator(N,rmin,rmax,ir)
        {
            theCache<double>().Register(this); //Can this move to the evaluator level?
        };
        EOrbital_HF_IBS(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
        : IrrepBasisSetImp<Evaluator>(ir)
        , Evaluator(es,ir,ltrim)
        {
            theCache<double>().Register(this);
        };


        virtual rFIT_CD_ABS* CreateCDFitBasisSet(const Structure*) const
        {
            return new EFit_IBS(Evaluator::Rescale(2.0, this->GetSymt()));   // the Coulomb fit: no mesh
        }
        virtual rFIT_SF_ABS* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
        {
            return new EFit_IBS(Evaluator::Rescale(2.0/3.0, this->GetSymt()), *cl, mp);
        }
        //! This block's Coulomb-fit exponents (its own x 2), for \c CoulombFitLadder's merged rungs.
        rvec_t CoulombFitExponents() const {return 2.0*this->es;}
        //! A Coulomb fit of THIS block's symmetry over \a es (verbatim): the merged rungs, on the first (s) block.
        rFIT_CD_ABS* MakeCoulombFitOver(const rvec_t& es) const {return new EFit_IBS(Evaluator(es, this->GetSymt()));}
        //! A power-m Coulomb fit (\a es at half-powers \a ms) of THIS block's symmetry, on the first (s) block: the family's own rungs.
        rFIT_CD_ABS* MakeCoulombFitOver(const rvec_t& es, const std::vector<int>& ms) const {return new EFit_IBS(Evaluator(es, ms, this->GetSymt()));}
        int BlockL() const {return Evaluator::Getl();}

        virtual std::ostream&  Write(std::ostream& os) const
        {
            Orbital_1E_IBS<Evaluator>::Write(os);
            Evaluator::Write(os);
            return os << std::endl;
        }

    };

    BasisSet_HF(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        for (auto ir:ec.GetIrreps())
            { auto* b=new EOrbital_HF_IBS(N,remin,remax,ir); Insert(b); itsBlocks.push_back(b); }
     
    }
    BasisSet_HF(const rvec_t& es, const ElectronConfiguration& ec, size_t ltrim=0)
    {
        for (auto ir:ec.GetIrreps())
            { auto* b=new EOrbital_HF_IBS(es,ir,ltrim); Insert(b); itsBlocks.push_back(b); }

    }
    // INDEPENDENT per-l exponent lists: each occupied irrep gets its OWN list \c esByL[l], verbatim (no
    // ExponentScaler auto-trim -- the trim is wrong for pseudopotential valence work, where the tight core
    // exponents it assumes are absent).  The caller must supply exponents for every occupied l.
    BasisSet_HF(const std::map<int,rvec_t>& esByL, const ElectronConfiguration& ec)
    {
        for (auto ir:ec.GetIrreps())
        {
            const int l=Symmetry::Atom::Getl(ir);
            auto it=esByL.find(l);
            // Exponents are REQUIRED for every occupied angular momentum.  A missing l is a caller error
            // (an incomplete recipe), not an internal invariant -- and it is reachable from untrusted CLI
            // input (CLIapps/valgen), so THROW rather than only assert: an unguarded end() deref here builds
            // an IBS from a garbage exponent list and aborts with std::bad_alloc far from the real cause.
            assert(it!=esByL.end());
            if (it==esByL.end())
                throw std::runtime_error("Atom BasisSet_HF: the electron configuration occupies angular "
                    "momentum l=" + std::to_string(l) + " but the basis recipe supplied no exponents for it. "
                    "Add a shell for l=" + std::to_string(l) + ".");
            { auto* b=new EOrbital_HF_IBS(it->second,ir,0); Insert(b); itsBlocks.push_back(b); } // ltrim=0: this l's full list
        }
    }

    virtual rFIT_CD_ABS* MakeCoulombFitRung(int rung) const override
    {
        if (rung<0 || rung>=kCoulombFitRungs)
            throw std::out_of_range("Radial::BasisSet_HF::MakeCoulombFitRung: rung "+std::to_string(rung)+" outside [0, "+std::to_string(kCoulombFitRungs)+")");
        if (itsBlocks.empty()) throw std::logic_error("Radial::BasisSet_HF::MakeCoulombFitRung: a basis with no blocks has no Coulomb fit");
        if (rung==0) return itsBlocks.front()->CreateCDFitBasisSet(nullptr);   // the hub's own fit, verbatim
        std::vector<rvec_t> perBlock;
        for (const auto* b : itsBlocks) perBlock.push_back(b->CoulombFitExponents());
        // A family whose fit functions carry their own radial power (Slater) fits each l-shell with its own shape, r^{2l} e^{-2 zeta r};
        // the others (Gaussian: its rung 0 already meets eps on Xe and U, 2026-11-14) merge every block's exponents into one s set.
        if constexpr (std::is_constructible_v<Evaluator, const rvec_t&, const std::vector<int>&, const sym_t&>)
        {
            std::vector<int> ls;
            for (const auto* b : itsBlocks) ls.push_back(b->BlockL());
            const auto [es, ms]=ShellPowerCoulombFit(perBlock, ls, rung);
            return itsBlocks.front()->MakeCoulombFitOver(es, ms);
        }
        else
            return itsBlocks.front()->MakeCoulombFitOver(MergedCoulombFitExponents(perBlock, rung));
    }
private:
    std::vector<const EOrbital_HF_IBS*> itsBlocks;   //!< the blocks, as built (owned by BasisSetImp): the Coulomb-fit ladder reads their exponents
};

// NR: 1E and HF, but no DFT
template <is1E_HF_Evaluator Evaluator> class BasisSet_1E_HF
    : public virtual Real_BS
    , public BasisSetImp<double>
{
public:
    // Define the Orbital IBS using mixins.
    class Orbital_1E_HF_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_ERI4_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
    {
    public:
        Orbital_1E_HF_IBS(size_t N, double rmin, double rmax, const sym_t& yl)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(N,rmin,rmax,yl)
        {
            theCache<double>().Register(this);
        };
        Orbital_1E_HF_IBS(const rvec_t& es, const sym_t& yl)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(es,yl)
        {
            theCache<double>().Register(this);
        };
         virtual std::ostream&  Write(std::ostream& os) const
        {
            Orbital_1E_IBS<Evaluator>::Write(os);
            Evaluator::Write(os);
            return os << std::endl;
        }
    };

    BasisSet_1E_HF(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new Orbital_1E_HF_IBS(N,remin,remax,ir));  
     
    }
    BasisSet_1E_HF(const rvec_t& es, const ElectronConfiguration& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new Orbital_1E_HF_IBS(es,ir));  
     
    }

};

template <is1E_HF_Evaluator LEvaluator, is1E_Evaluator SEvaluator> class BasisSet_RKB
    : public virtual Real_BS
    , public BasisSetImp<double>
{
public:
    class EOrbital_RKBL_IBS 
        : public Orbital_RKBL_IBS<LEvaluator>
        , public Orbital_ERI4_IBS<LEvaluator>
        , public IrrepBasisSetImp<LEvaluator>
        , public LEvaluator
    {
    public:
        EOrbital_RKBL_IBS(size_t N, double emin, double emax, const sym_t& irrep)
        : IrrepBasisSetImp<LEvaluator>(irrep)
        , LEvaluator(N,emin,emax,irrep)
        {
            theCache<double>().Register(this); //Can this move to the evaluator level?
        };
        EOrbital_RKBL_IBS(const rvec_t& es, const sym_t& irrep, size_t ltrim=0)
        : IrrepBasisSetImp<LEvaluator>(irrep)
        , LEvaluator(es,irrep,ltrim)
        {
            theCache<double>().Register(this);
        };

        virtual std::ostream&  Write(std::ostream& os) const
        {
            Orbital_RKBL_IBS<LEvaluator>::Write(os);
            LEvaluator::Write(os);
            return os << std::endl;
        }
    };
    class EOrbital_RKBS_IBS 
        : public Orbital_RKBS_IBS<SEvaluator>
        , public IrrepBasisSetImp<SEvaluator>
        , public SEvaluator
    {
    public:
        EOrbital_RKBS_IBS(size_t N, double emin, double emax, const sym_t& irrep)
        : IrrepBasisSetImp<SEvaluator>(irrep)
        , SEvaluator(N,emin,emax,irrep) //fix κ=-1, l=0
        {
            theCache<double>().Register(this); //Can this move to the evaluator level?
        };
        EOrbital_RKBS_IBS(const rvec_t& es, const sym_t& irrep, size_t ltrim=0)
        : IrrepBasisSetImp<SEvaluator>(irrep)
        , SEvaluator(es,irrep,ltrim)
        {
            theCache<double>().Register(this);
        };

        virtual std::ostream&  Write(std::ostream& os) const
        {
            Orbital_RKBS_IBS<SEvaluator>::Write(os);
            SEvaluator::Write(os);
            return os << std::endl;

        }
    };
    class EOrbital_RKB_IBS 
        : public virtual Orbital_RKB_IBS<double>
        , public virtual ::qchem::BasisSet::Orbital_ERI4_IBS<double>
        , public Orbital_RKB_HF_IBS_Imp<double>
        , public ::qchem::BasisSet::IrrepBasisSetImp<double>
        , public Symmetry::IdentityRepresentation<double>   // radial κ block: identity under ops fixing the nucleus (S3)
    {
    public:
        EOrbital_RKB_IBS(size_t N, double remin, double remax, const sym_t& irrep)
        : Orbital_RKB_HF_IBS_Imp(
                    new EOrbital_RKBL_IBS(N,remin,remax,irrep),
                    new EOrbital_RKBS_IBS(N,remin,remax,irrep)
                )
        , IrrepBasisSetImp<double>(irrep)
        , Symmetry::IdentityRepresentation<double>(Symmetry::Stabiliser::FixingPoint(rvec3_t(0,0,0)))
        {};
        EOrbital_RKB_IBS(const rvec_t& es, const sym_t& irrep, size_t ltrim=0)
        : Orbital_RKB_HF_IBS_Imp(
                    new EOrbital_RKBL_IBS(es,irrep,ltrim),
                    new EOrbital_RKBS_IBS(es,irrep,ltrim)
                )
        , IrrepBasisSetImp<double>(irrep)
        , Symmetry::IdentityRepresentation<double>(Symmetry::Stabiliser::FixingPoint(rvec3_t(0,0,0)))
        {};

        virtual size_t GetNumFunctions() const {return Orbital_RKB_IBS_Imp<double>::GetNumFunctions();}
    };
    BasisSet_RKB(size_t N, double ermin, double ermax, const ElectronConfiguration& ec)
    {
        for (auto irrep:ec.GetIrreps())
            Insert(new EOrbital_RKB_IBS(N,ermin,ermax,irrep));

    }
    BasisSet_RKB(const rvec_t& es, const ElectronConfiguration& ec, size_t ltrim=0)
    {
        for (auto irrep:ec.GetIrreps())
            Insert(new EOrbital_RKB_IBS(es,irrep,ltrim));
    }

};


}} //namespaces 

