class ElectronConfigurationTests;

class DiracIntegralTests;
class SlaterRadialIntegralTests;
class Cache4Tests;

class RunDiagnosticsTests;   // src/Calculation/tests/RunDiagnostics.C -- the T2/T3/T4
                             // detector rules, tested on SYNTHETIC trajectories rather than by
                             // coaxing a full SCF into the outcome under test.

class TransitionDensityTests;   // src/Response/tests -- the finite-difference response-kernel ORACLE reads the
                                // AO blocks of a transition density (doc/LinearResponsePlan.md ruling D6);
                                // the TransitionDensity face never hands its matrices out.
class DeckTests;                // src/Calculation/tests/Deck.C -- an inline-structure deck and its named twin must build the
                                // SAME restart fingerprint (SolidCalculation::Fingerprint is private).
class ResponseFacadeTests;      // src/Response/tests -- the same oracle, driven through a converged
                                // qchem::MoleculeCalculation (its Hamiltonian / wave function / basis).
class ForceFacadeTests;         // src/Calculation/tests/ForceFacade.C -- a term's energy at moved atoms with D frozen (the per-term
                                // finite difference of doc/GeometryOptimizationPlan.md §6) needs MoleculeCalculation's Hamiltonian.
class PointJetsTests;           // src/BasisSet/Gaussian/tests/M_GeometryDerivatives.C -- fills a PointJets with analytic functions to check its
                                // jets, ρ's jet and the density pullback against finite differences, before any family's JetsOn exists
                                // (GeometryOptimizationPlan §13.2: JetsOn's gradient source awaits the user's ruling).
