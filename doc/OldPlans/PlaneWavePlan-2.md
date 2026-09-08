# Code Review 
Following the pipline Si SCF DFT test, line 1016 in UnitTests/PlaneWaveDFTUT.C




We need to copy src/BasisSet/Molecule/Evaluators/Evaluator.C over to src/BasisSet/lattice_3D/Evaluators/Evaluator.C and adjust the namespace and module name accordingly. Then we are free to (grudginly) tune these interfaces to meet the requirements for all the Lattie_3D IBS classes.

In src/BasisSet/Molecule/IrrepBasisSet.C we can see that the Orbital_1E_IBS class supports both single element evaluators (is1E_Evaluator)  and full matrix Evaluators (isM_1E_Evaluator) for Make{Overlap,MakeKinetic, MakeNuclear}.  But the Molecule example is either all single element or all matrix.  For PW it appears we have a mixture.  Can we function foward inline single element calls for all three Make{Overlap,MakeKinetic, MakeNuclear}?  We can still capture the diagonal flavour for Overap and Kinetic in Orbital_1E_IBS<is1E_Evaluator>. SOmething like:


At the end of each evaluator class we should have a statement like:
```c++
static_assert(is1E_Evalualtor<Evaluator>);
```
If the is1E_Evalualtor is impossible and we must use isM_1E_Evalualtor that should also be fine.  The goal is put the least possible burden on Evaluator developer, which means is1E_Evalualtor is prefered.  Then our collection of IBS basis sets is maximally extensible.

