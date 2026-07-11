// File: SCFParams.C Parameters used for controling SCF iteration and convergence.
export module qchem.SCFParams;
import qchem.Types;

namespace qchem {

// Defaults are the proven molecular-HF recipe (matches M_HF_U / the bridge): they converge a
// neutral molecule (e.g. water/dzvp) out of the box, so a facade caller can write
// `calc.Converge({.NMaxIter=60})` and override only what they care about.  Positional aggregate
// init (`{20,1e-4,...}`) and per-field assignment (`par.MinΔρ=...`) both still work unchanged --
// default member initializers do not disqualify the aggregate.
export struct SCFParams
{
    size_t NMaxIter        = 20;     //Max allowed number of iterations
    double MinΔρ           = 1e-4;   //Minimum delta in charge density for convergence.
    double MinΔFD          = 1e-7;   //Minimum delta in [F,D] (|[F,D]-[F,D]_old|) for convergence.
    double MinΔE           = 1e30;   //Minimum RELATIVE total-energy change |ΔE/E| for convergence (default
                                     //  off).  The physical gate for a NON-variational SCF (density-fit /
                                     //  GPW collocation) whose [F,D]/Δρ limit-cycle above the fit floor while
                                     //  the total energy is settled -- see doc/GPWPlan.md.
    double MinVirial       = 1e-13;  //Minimum error Virial ratio -V/K.  (1e-13 => effectively off; the
                                     //  textbook -V/K=2 virial is not gated for molecules.)
    double MinFD           = 1e-5;   //Minimum error from SCF accelerator.  i.e. [F,D] (orbital gradient).
    double StartingRelaxRo = 1.0;    //relaxation for mixing Ro.  Dynamically adjusted during iterations.
    double MergeTol        = 1e-4;   //Merge eigen levels (like Px,Py Pz) that are equal within +/- MergeTol
    bool   Verbose         = false;  //Display iteration details.
};



} // namespace qchem