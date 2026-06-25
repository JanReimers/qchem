# QCMesh project-module (library) upgrade plan

- 3D meshes are built from Radial(x)Angular parts with numerous options for each component.  Yet another
use of abstractions to handle the combinitorics.
- RadialMesh is 1D and has its own interface.
- Angular mesh is 3D and uses the Mesh interface. 
- So within the qcMesh project module you might ask: how do you comibine a RadialMesh and AngularMesh and get a ready to use
  instance of a Mesh?  Answer: You can't!!!!  That functionality go refactored over into qcStructure which knows how to make
  3D meshes for Molecules using Varonoi polehedra.
- Is this a design flaw?  Maybe. 
- SO now we have a couple of TODO items from the PlaneWave project:
1. **qcMesh owns quadrature**: PlaneWave_IBS's hand-rolled `Integral(f)=Sum Omega/Npts*f` + UniformGrid
  belongs in qcMesh (Mesh carries points+weights; integral = Sum w_i f(r_i)). Then migrate all PW IBSs.
2. **Reconcile mesh ownership**: NEW model = the basis owns the mesh (correct); OLD code has fitted
  potentials owning the Mesh (likely wrong) — revisit after the qcMesh refactor.
  If we decide to refactor this for CD and Vxc fitting then this might be a big job.  But for atoms at least we have a DFT oracle now so the risk is low.
3. High level user (UnitTests/PlaneWaveDFTUT.C lines 1055-1062) needs a way to specify the real space unit cell integrationg grid.
- There work in qcMesh and qcStructure and bleed into other project-modules.  WHich is fine.
- In addition while we are in qcMesh I would like to a couple of other things:
1. The Mesh interface exposes iterator for RW={rvec3_t r,double weight} pairs.  They are also stored as pairs and this is precicesly backwards from what we want because the R and W parts are used in totally different places.  Rs are used in op(r) evaluators, and Ws are used in the integrations algos.  They need separat storage and separate iterators.
2. We should consider harmonizing with the iterator pattern used for Structure.
3. The Clone() function can probably be removed.
4. The complex MeshParams struct shoudl probably be replaced with json argument to the factory.  That way the user not required
to define paramaters that are not even used.
5. The molecular mesh in src/Structure/Internal/MoleculeMesh.C has a bug somewhere.  If you create a dimer with the two atoms
right on top of each other, the integrals should evaluate to the same values as integrating an atom mesh.  They don't.
6. For FLAPW we will need a unit cell analogue of the molecular mesh.  Varonoi polyhedra chopped off at the cell walls.
7. The Atom BSpline evaluator has its own little GL integrator in src/BasisSet/Atom/Evaluators/BSpline/Internal/GLQuadrature.C using the fortran code from Numerical recipes: src/BasisSet/Atom/Evaluators/BSpline/Internal/gauleg.f Can we convert the fortran
to c++ with blaze vectors and incorporate into qcMesh.  Either way the BSpline and PW basis should use the same 1D GL code if possible.
8. qcMesh has ScalarFunction and VectorFunction ... these are very useful interfaces.  BUT I made the mistake of including virtual mat_t<T> operator() (const Mesh&   ) const  ; in both of them ... this everything that wants to inherit and interface for virtual vec_t<T> operator() (const rvec3_t&) const=0; gets Mesh stuff dragged in there with it.  ISP violation.  If you touch Mesh.C the project rebuilds ... un SOLID.  We need to seprate and scalar an vector operator() (const rvec3_t&) const=0 interfaces separated from Mesh
9. src/Mesh/MeshIntegrator.C is over engineered ... all the class does is hold a Mesh* which by itself it useless.  These should just be plain old functions that take const Mesh& as the first paramater.
10. We should remove all the repulsion integrals from src/Mesh/MeshIntegrator.C ... they can't work.  I didn't read Hartree!

OK that is it for now ... I usually think of more things once the work starts!!  We are doing agile not waterfall :)