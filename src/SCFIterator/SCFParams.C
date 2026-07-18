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
    double KerkerG0        = 0.0;    //Kerker density-mixing screening wavevector G0 (a.u.^-1).  0 (default) =
                                     //  OFF -> classic linear density-matrix mixing (atoms/molecules unchanged).
                                     //  >0 -> preconditioned rho-mixing on the periodic (GPW/PW) path: damp the
                                     //  low-G charge-transfer slosh by G^2/(G^2+G0^2) (breaks the ionic-crystal
                                     //  limit cycle -- see doc/GPWPlan.md).  Typical G0 ~ 1.0.  dcmplx path only.
    bool   UseMOM          = false;  //Maximum Overlap Method: occupy the orbitals with the largest overlap onto a
                                     //  FIXED reference occupied subspace instead of the lowest eigenvalues.  OFF
                                     //  (default) -> plain aufbau (atoms/molecules unchanged).  ON -> occupied-
                                     //  subspace continuity, the fix for a diving diffuse virtual being aufbau-
                                     //  captured (the NaF Γ occupation-swap instability -- see doc/GPWPlan §0b″).
    int    PulayDepth      = 0;      //Density-face Pulay (density-DIIS) history depth on the periodic ρ̃ path.  0
                                     //  (default) = plain Kerker/linear mixing.  >0 = Kerker-preconditioned Pulay
                                     //  mixing keeping this many (ρ̃_in, ρ̃_out) pairs (VASP/QE/CP2K scheme; the
                                     //  quasi-Newton cure for the charge-transfer slosh -- see doc/SCFStrategyPlan.md).
    int    PulayStart      = 0;      //Prime with plain Kerker for this many iterations before engaging Pulay (only
                                     //  used when PulayDepth>0).  History-based mixing is UNSTABLE far from the fixed
                                     //  point (residuals not yet in the linear-response regime), so descend first,
                                     //  then accelerate -- the density-side of the ladder hand-off.  0 = immediate.
    int    MOMStartIter    = 10;     //Delayed-IMOM reference-capture iteration (only used when UseMOM).  Run plain
                                     //  aufbau for this many fills so the SCF descends to the physical fixed point,
                                     //  THEN capture the occupied subspace ONCE and hold it fixed (the seed is
                                     //  mid-transient, so capturing at iteration 0 anchors garbage).
};



} // namespace qchem