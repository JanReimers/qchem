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
    //! \brief What a solid run's accelerator chain is made of, defaulted to the GPW production recipe.
    //!
    //! COMPOSES the public \c DIISParams / \c GDMParams / \c ScheduleSignal rather than restating their
    //! fields.  An earlier version flattened them into renamed scalars, which was a drift trap: a new knob
    //! on \c DIISParams would silently not reach this door.  If you are tuning DIIS, you are setting
    //! \c DIISParams -- there should be exactly one spelling of that.
    struct SolidAcceleratorOptions
    {
        DIISParams diis {8, 0.1, 1e-10, 1e-9};   //!< Nproj, then the [F,D] gates (NOT energy -- see DIISParams).
        GDMParams  gdm  {1.0, 0.1};              //!< standalone GDM.
        //! Ladder only: the GDM rung takes a larger trust radius than standalone, because it is polishing a
        //! tail DIIS has already brought close rather than descending from far away.
        GDMParams      ladderGdm      {1.0, 0.5};
        double         ladderEThresh  = 1e-8;    //!< STALL trigger: energy still moving by more than this.
        int            ladderStall    = 5;       //!< ...for this many non-improving steps.
        double         ladderFloor    = 1e-8;    //!< low backstop on the stall trigger.
        double         ladderSwitchAt = 1e-6;    //!< TAIL trigger: hand off to GDM once the signal is below this.
        //! \c EnergyChange, not the molecular \c Error default: on a solid, [F,D] is contaminated by
        //! charge-transfer slosh and may never settle while the total energy does (see \c ScheduleSignal).
        ScheduleSignal ladderSignal   = ScheduleSignal::EnergyChange;
    };

    //! \brief Build the complex (Bloch) accelerator a solid run asks for.  The \c cSCFAccelerator twin of
    //! the factory above, which did not exist: \c cSCFAccelerator is a public alias but every concrete
    //! complex accelerator is \c .Internal., so a solid run could only be driven from a unit test.
    //! \c Ladder is the ionic-crystal recipe -- DIIS until it runs out of steam, then GDM, switching on
    //! ENERGY CHANGE rather than the residual.
    cSCFAccelerator* Factory(Type, const SolidAcceleratorOptions& = {});
}

