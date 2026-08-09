// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.SCFAccelerator.Factory;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{
    //! \c Null = no acceleration (plain damped iteration) -- a real policy for a run whose mixer does the
    //! work (NaF's pure damped Kerker), and the honest way to A/B an accelerator against none.
    enum class Type {DIIS,GDM,Ladder,Null};

    //! The molecular (\c double) front door.  Stringly-typed by history: the json keys came from test
    //! scaffolding.  The typed \c SolidAcceleratorOptions below is the pattern to prefer for new callers.
    SCFAccelerator* Factory(Type,const nlohmann::json& js);

    //=== The SOLID (dcmplx / Bloch) front door ========================================================
    //! \brief Typed knobs for a complex accelerator chain, defaulted to the GPW production recipe.
    //!
    //! TYPED, not json: this door exists for a FACADE (\c SolidCalculation), and the concrete
    //! \c DIISParams / \c GDMParams / \c ScheduleSignal live in \c .Internal. modules that no facade may
    //! import across a library boundary.  Naming the few knobs here keeps the internals internal -- the same
    //! reason \c AcceleratorOptions exists one layer up for the molecular path.
    //!
    //! NB \c fdMax / \c fdMin gate on the RESIDUAL \f$[F,D]\f$, not the energy (see \c DIISParams: the old
    //! "EMax" spelling was a documented trap).
    struct SolidAcceleratorOptions
    {
        size_t nProj  = 8;         //!< DIIS history depth.
        double fdMax  = 0.1;       //!< DIIS starts extrapolating once [F,D] < fdMax.
        double fdMin  = 1e-10;     //!< DIIS stops once [F,D] < fdMin.
        double svTol  = 1e-9;      //!< DIIS bails when min singular value of B < svTol.
        double gdmFDMax = 1.0;     //!< GDM engages once [F,D] < this.
        double gdmTrust = 0.1;     //!< GDM trust radius (radians per step).
        // --- Ladder only: when to hand DIIS over to GDM ---
        double ladderGdmTrust = 0.5;
        double ladderEThresh  = 1e-8;
        int    ladderStall    = 5;
        double ladderFloor    = 1e-8;
        double ladderSwitchAt = 1e-6;
    };

    //! \brief Build the complex (Bloch) accelerator a solid run asks for.  The \c cSCFAccelerator twin of
    //! the factory above, which did not exist: \c cSCFAccelerator is a public alias but every concrete
    //! complex accelerator is \c .Internal., so a solid run could only be driven from a unit test.
    //! \c Ladder is the ionic-crystal recipe -- DIIS until it runs out of steam, then GDM, switching on
    //! ENERGY CHANGE rather than the residual.
    cSCFAccelerator* Factory(Type, const SolidAcceleratorOptions& = {});
}

