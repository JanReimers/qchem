// File: ElectronConfigurations/ElectronConfiguration.C Interface for and electron configuration.
module;
#include <set>
export module qchem.ElectronConfiguration;
export import qchem.Symmetry.Irrep;

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
};