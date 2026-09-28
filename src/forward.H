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
class ResponseFacadeTests;      // src/Response/tests -- the same oracle, driven through a converged
                                // qchem::Calculation (its Hamiltonian / wave function / basis).
