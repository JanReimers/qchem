// File: ElectronConfigurations/ElectronConfiguration.C Interface for and electron configuration.
module;
#include <set>
export module qchem.ElectronConfiguration;
export import qchem.Symmetry.Irrep;

namespace qchem {

export class ElectronConfiguration
{
public:
    // Define how symmetries are to be ordered
    static constexpr auto cmp = [](sym_t a, sym_t b) 
    {
            return a->SequenceIndex() < b->SequenceIndex();
    }; 
    using syms_t=std::set<sym_t,decltype(cmp)>;

    virtual ~ElectronConfiguration() {};
    virtual int    GetN(const Irrep&) const=0;
    virtual syms_t GetIrreps() const=0;
    virtual void   Display() const=0;
    // When true, the wave function fills the globally-lowest orbitals across all irreps each
    // iteration (a molecular aufbau) and GetN gives the TOTAL per spin channel; when false the
    // per-irrep GetN count is taken as fixed (atoms / hand-set occupations).
    virtual bool   UsesAufbau() const {return false;}
    // When true (a metal), the Bloch k-blocks share ONE chemical potential across the BZ mesh: the fill
    // solves a single μ on Σ_k w_k Σ_i g_i f_i = GetN (the whole-mesh total per spin channel) and charge
    // sloshes between k-points, instead of pinning each block to a fixed per-block count.  Mutually
    // exclusive with UsesAufbau (a Bloch block IS an irrep; there is no cross-irrep aufbau to run).
    // doc/GPWPlan1.md item 3.
    virtual bool   UsesGlobalFermi() const {return false;}
};
} // namespace qchem