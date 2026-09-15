// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.

/*!
\page basisset_taxonomy The basis-set taxonomy: one road map for every basis type and every engine

\tableofcontents

This code holds, in ONE electronic-structure framework, radial bases for atoms (Slater, Gaussian,
B-spline, non-relativistic and Dirac), Cartesian / spherical Gaussians for molecules, plane waves and
Gaussian-and-plane-waves (GPW) for crystals, delta (grid) bases for density and potential fits, and has
slots for augmented plane waves, numerical atomic orbitals and more.  As far as we know, nobody has tried
to combine that many basis families and integral engines under one set of interfaces before, so the
organising principle is written down here in one place, with the definitions a chemist needs to read the
placement table at the end.  The design ruling itself is `doc/BasisSetTaxonomyPlan.md` (a RECORD file:
executed 2026-09-13); this page is its user-facing statement and is meant to stay current.

If you know what an irrep is (you have used a character table to say that water's O 2p\f$_z\f$ orbital is
\f$a_1\f$ and its two H 1s orbitals combine into \f$a_1\oplus b_2\f$), you have everything needed.  Every
other term is defined where it first appears and collected again in the \ref bst_glossary "glossary".

\section bst_one_idea 1. The one idea: a basis block is a carrier space for one irrep of G

Let \f$G\f$ be the **symmetry group of the Hamiltonian**: the set of all spatial (and, when relevant, spin)
operations \f$g\f$ that leave \f$\hat H\f$ unchanged, \f$g\hat H g^{-1}=\hat H\f$.  For a free atom that is
\f$O(3)\f$; for water it is \f$C_{2v}\f$; for a crystal it is (at least) the group of lattice translations
\f$T\f$.  Because \f$\hat H\f$ commutes with every \f$g\in G\f$, \f$\hat H\f$ is **block-diagonal** in a basis
whose functions are sorted by irrep: a matrix element \f$\langle\phi|\hat H|\psi\rangle\f$ vanishes unless
\f$\phi\f$ and \f$\psi\f$ transform as the same irrep \f$\Gamma\f$ (and the same row of it).  Written out,

\f[
   \hat H \;=\; \bigoplus_{\Gamma\in\hat G}\; \hat H_\Gamma ,
   \qquad\qquad
   \mathrm{BasisSet} \;=\; \bigoplus_{\Gamma}\; \mathrm{IrrepBasisSet}_\Gamma .
\f]

The symbol \f$\oplus\f$ is the **direct sum**: the whole space is the concatenation of independent
sub-spaces, one per irrep, with no matrix elements between them.  In the code that sentence is literal:
`tBasisSet<T>` is a list of `IrrepBasisSet<T>` blocks, and the SCF solves one block at a time.

Each block is a **carrier space** for its irrep.  A representation of \f$G\f$ is a rule assigning a matrix
\f$D(g)\f$ to every group element; the *carrier space* is the vector space those matrices act on -- the set
of functions that get mixed among themselves, and only among themselves, by the operations of \f$G\f$.  An
`IrrepBasisSet` is therefore not "a basis with symmetry \f$\Gamma\f$" in the loose sense; it is precisely a
set of functions spanning copies of the carrier space of \f$\Gamma\f$.

Two consequences dissolve most of the historical confusion about what belongs where:

-# **Carriers are not invariants.**  \f$Y_{lm}\f$ is not spherically symmetric; \f$e^{i\mathbf k\cdot\mathbf r}\f$
   is not translation-invariant; a SALC is not point-group-invariant.  Each *transforms as* an irrep.
   Invariance was never the criterion for membership in a block -- transforming correctly is.
-# **\f$\mathbf k\f$ is an irrep label of the translation group \f$T\f$**, exactly as \f$l\f$ labels the
   irreps of \f$O(3)\f$ and \f$a_1\f$ labels an irrep of \f$C_{2v}\f$.  \f$T\f$ is abelian, so all its irreps
   are one-dimensional; the irrep labelled \f$\mathbf k\f$ sends the translation by \f$\mathbf R\f$ to the
   phase \f$e^{-i\mathbf k\cdot\mathbf R}\f$, and its carriers are the Bloch functions.  Bloch's theorem is
   nothing but the statement that the projection operator onto that irrep,
   \f[
      \hat P_{\mathbf k} \;=\; \frac{1}{N}\sum_{\mathbf R} e^{-i\mathbf k\cdot\mathbf R}\,\hat T_{\mathbf R},
   \f]
   block-diagonalises \f$\hat H\f$.  So \f$\mathbf k\f$ is neither "in the basis" nor "in the eigenfunctions":
   it is the block label, and the eigenfunctions inherit it because \f$\hat H\f$ is block-diagonal in it.
   The code says so: `tGPW_IBS(const UnitCell&, const sym_t& irrep, ...)` -- the Bloch symmetry IS the
   \f$\mathbf k\f$-label.

\section bst_axes 2. Three orthogonal axes and a role

Every basis in the table below is placed by answering three independent questions, plus one about what
it is used for.  The axes are orthogonal: fixing one says nothing about the others.

| axis | question | values today | where it lives in the code |
|---|---|---|---|
| **G** | which group's irreps label the blocks? | \f$1\f$, \f$P\f$, \f$O(3)\f$, \f$O(3)^*\f$, \f$T\f$; later \f$T\rtimes P\f$, \f$P^*\f$, \f$T^d\f$ | `qcSymmetry/{Atom, Molecule, Lattice_3D}` -- laid out by G |
| **Family** | which elementary functions does the engine integrate? | Gaussian (primitive / contracted, Cartesian / spherical); Slater; B-spline; exponential (plane wave); delta; composite (APW) | the engine libraries `qcRadial_BS`, `qcGaussian_BS`, `qcPlaneWave_BS` |
| **Construction** | how does the seed become an irrep carrier? | subduction (A) or induction (B) -- derived from (G, family), \ref bst_construction "section 3" | `qcSymmetry` (SALC, Bloch phases, Fold) |
| **Role** | what is the basis for? | Orbital (1E, HF, DFT, DHF, PP faces); Fit (charge density, scalar field) | the `qcBasisSet` core faces |

\subsection bst_groups 2.1 The groups on the G axis

- \f$\mathbf 1\f$ -- the **trivial group**, one element, one irrep.  "No symmetry" is a row in the table,
  not a separate category: a non-symmetric molecule simply has one block.
- \f$P\f$ -- a **point group** (\f$C_{2v}\f$, \f$D_{3h}\f$, \f$O_h\f$, ...): rotations, reflections and
  improper rotations that fix a point.  The chemist's home ground.
- \f$O(3)\f$ -- the **full rotation group of three-dimensional space** (all rotations about a point, proper
  and improper, i.e. including inversion): the symmetry of a free atom.  Its irreps are labelled by the
  angular momentum \f$l=0,1,2,\dots\f$ (and parity), have dimension \f$2l+1\f$, and their carrier space is
  spanned by the spherical harmonics \f$Y_{lm}\f$, \f$m=-l\dots l\f$.  Note that \f$m\f$ labels the ROW
  within the irrep, not the irrep: for an \f$O(3)\f$-symmetric \f$\hat H\f$ Schur's lemma gives
  \f$\langle nlm|\hat H|n'l'm'\rangle=\delta_{ll'}\delta_{mm'}\,h^{(l)}_{nn'}\f$ with \f$h^{(l)}\f$ independent of
  \f$m\f$, so the \f$2l+1\f$ rows are identical blocks and the atom solves ONE radial problem per \f$l\f$
  (the \f$Y_{lm}\f$ "is carried by the irrep, not the function" -- `ImplicitAngular_IBS`).  This is a
  statement about \f$G\f$ being \f$O(3)\f$, not about atoms: an open-shell atom whose unpaired electrons are
  placed in specific \f$m\f$'s (the maximally stretched determinant `Atom_EC` builds -- carbon \f$p^2\f$ as
  \f$m\in\{-1,0\}\f$ unpaired, \f$\{+1\}\f$ paired) has a Fock operator that is only AXIALLY symmetric, so
  \f$G\f$ is the smaller group, \f$m\f$ (there, the unpaired and paired \f$m\f$-sets, `YFactory(l, ms)`)
  becomes an irrep label, and those blocks legitimately differ.  Every finite point group is a subgroup
  of \f$O(3)\f$, which is why crystal-field splitting is a *branching rule* (below).
- \f$O(3)^*\f$, \f$P^*\f$, \f$(T\rtimes P)^*\f$ -- the **double groups**.  A spin-\f$\tfrac12\f$ particle picks up a
  sign under a \f$2\pi\f$ rotation, so a rotation by \f$2\pi\f$ is no longer the identity; the group is
  doubled and gains extra ("spinor") irreps.  For the atom these are the Dirac \f$\Omega_{\kappa m_j}\f$
  two-spinors that label the DHF blocks.
- \f$T\f$ -- the **translation group** of a lattice: \f$\{\hat T_{\mathbf R}\}\f$ for all lattice vectors
  \f$\mathbf R\f$.  Abelian; irreps labelled by \f$\mathbf k\f$ in the Brillouin zone; carriers = Bloch
  functions.  \f$T^d\f$ with \f$d=1,2\f$ is the same thing for a polymer or a slab.
- \f$T\rtimes P\f$ -- a **space group**.  The symbol \f$\rtimes\f$ is the **semidirect product**: the group
  whose elements are pairs "rotation \f$R\in P\f$ followed by translation \f$\mathbf t\f$", written
  \f$\{R|\mathbf t\}\f$, with \f$T\f$ a *normal* subgroup on which \f$P\f$ acts by rotating the translation
  vectors, \f$\{R|\mathbf 0\}\{E|\mathbf t\}\{R|\mathbf 0\}^{-1}=\{E|R\mathbf t\}\f$.  It is not the direct
  product \f$T\times P\f$ because rotations and translations do not commute.  Its irreps are labelled by a
  star of \f$\mathbf k\f$-vectors plus an irrep of the little group of \f$\mathbf k\f$ -- the next G row the
  code will grow.
- \f$G_{\rm spatial}\times SU(2)\f$ -- the **direct product** with the spin rotation group.  \f$\times\f$ means
  the two factors commute and every irrep is a pair (spatial irrep, spin irrep).  This is the
  non-relativistic case, and it is *why* Pol/UnPol is a block structure at all: the two spin channels are
  the irreps of the \f$SU(2)\f$ factor.  `Symmetry::Irrep = (sym_t spatial) × (Spin ms)` encodes exactly
  this.  Spin-orbit coupling is the factor dissolving into \f$G\f$ -- the double-group rows above.

\subsection bst_family 2.2 The family axis: the analytic seed

The **seed** (or *analytic seed*) of a basis is the elementary function type the integral engine knows how
to integrate, before any symmetry adaptation has been done to it: a Gaussian \f$x^ay^bz^c e^{-\alpha r^2}\f$
on a centre, a Slater \f$r^{n-1}e^{-\zeta r}\f$, a B-spline segment, a plane wave \f$e^{i\mathbf q\cdot\mathbf r}\f$,
a delta function on a grid point.  A **family** is a seed type together with what can be built from it
(primitives, contractions, Cartesian or spherical shells...).  The family is the natural unit of an
integral **engine**: the McMurchie-Davidson recursions serve every Gaussian, the FFT serves every plane wave,
the radial quadrature serves every function of \f$r\f$ times a \f$Y_{lm}\f$.

Because the axes are orthogonal, **a family belongs to no G**.  `Radial/Evaluators/Gaussian` (atoms) and
`Gaussian/PG_Spherical` (molecules) are the SAME family with two engines; the atomic engine exists only
because \f$O(3)\f$ collapses every integral to a one-dimensional radial one.  "Gaussians are inefficient for
atoms through the molecular engine" is a statement about an engine, not about a family.

\subsection bst_role 2.3 The role axis

An **orbital** basis expands the one-electron states; its faces (`Orbital_1E_IBS`, `Orbital_HF_IBS`,
`Orbital_DFT_IBS`, `Orbital_DHF_IBS`, `Orbital_PP_IBS`) are the integrals a Hamiltonian term can ask it for.
A **fit** basis expands a *field* -- the charge density (`FIT_CD_ABS`, Coulomb metric) or a scalar potential
such as \f$v_{xc}\f$ (`FIT_SF_ABS`, overlap metric).  A quadrature IS a fit: a real-space integration grid is
a delta-function fit basis over its points, and "everything is a fit" = (integration grid) \f$\times\f$ (fit
family), two independent choices that are never hard-wired to each other (`doc/Pins.md`).

\section bst_construction 3. The construction axis is not free: subduce or induce

How does a seed (which usually carries some *other* symmetry) become a carrier of an irrep of the
Hamiltonian's \f$G\f$?  There are exactly two group-theoretical operations, and which one applies is fixed by
(G, family).  Both are standard representation theory; the chemist has met both without the names.

\subsection bst_subduction 3.1 Subduction (A): restricting a bigger group's carrier to G

If the seed is already adapted to a **bigger** group \f$G_{\rm big}\supseteq G\f$, its carrier space is a
representation of \f$G\f$ too -- just not necessarily an *irreducible* one.  **Subduction** (also called
*restriction*), written \f$G_{\rm big}\downarrow G\f$, is decomposing that representation into irreps of the
subgroup.  The recipe for doing so is a **branching rule**.

- Crystal-field splitting is the everyday example: a \f$d\f$ shell carries \f$l=2\f$ of \f$O(3)\f$; in an
  octahedral site, \f$O(3)\downarrow O_h\f$ gives \f$l{=}2\;\to\;t_{2g}\oplus e_g\f$.  Nothing new was built;
  the five functions were re-sorted into two blocks.
- In this code: a Gaussian shell on an atom carries \f$Y_{lm}\f$ (an \f$O(3)\f$ carrier); in a molecule it is
  subduced to the atom's **site group** \f$S_a\subseteq P\f$ (the point-group operations that leave that atom
  in place).  That is `Symmetry::Molecule::{Cartesian,Spherical}ShellRep`.
- Plane waves \f$e^{i\mathbf q\cdot\mathbf r}\f$ are carriers of the continuous translation group of all of
  \f$\mathbb R^3\f$; in a crystal they are subduced to the lattice translations, \f$\mathbb R^3\downarrow T\f$,
  which is precisely the \f$\mathbf q=\mathbf k+\mathbf G\f$ labelling: every wave with the same \f$\mathbf k\f$
  modulo a reciprocal-lattice vector \f$\mathbf G\f$ lands in the same block.

\subsection bst_induction 3.2 Induction (B): building G's carrier from a site-local seed

If the seed is **local** to a site with a small symmetry group \f$H\subset G\f$ (an AO on one atom, a delta on
one grid point), the operations of \f$G\f$ carry it to the other sites of its **orbit** (all the places
\f$G\f$ can move it to).  **Induction**, written \f$H\uparrow G\f$, is the representation of \f$G\f$ spanned by
the seed together with all its images; it has dimension \f$|G|/|H|\times\dim(\text{seed})\f$ and is
decomposed into irreps of \f$G\f$ by applying \f$G\f$'s projection operators.  (The formal theory is
Frobenius reciprocity; the practical recipe is the projection operator every group-theory course teaches.)

- The chemist's example: the two H 1s orbitals of water.  Each is local to its own H (site group
  \f$C_s\f$); \f$C_{2v}\f$ swaps them; inducing \f$C_s\uparrow C_{2v}\f$ over the orbit of the two H atoms
  gives \f$a_1\oplus b_2\f$ -- the symmetry-adapted linear combinations (SALCs) \f$1s_A\pm 1s_B\f$.
- In this code, molecules: AO shells on equivalent atoms \f$\to\f$ SALCs (site group \f$\uparrow P\f$),
  built in `Gaussian/Point/SymmetryAdaptedBasisSet`.
- Crystals: an AO on one cell atom \f$\to\f$ the **Bloch sum** \f$\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}
  \phi(\mathbf r-\mathbf R)\f$, which is \f$1\uparrow T\f$ evaluated with the projector of section 1.  This
  is what `Gaussian/Lattice/LatticeSum1E` and the GPW engine do.
- Deltas on a uniform grid \f$\to\f$ Bloch sums of deltas, which are plane waves restricted to the grid.
  **The uniform delta basis is the plane-wave basis's adjoint**: both carry \f$T\f$, one by induction, one
  by subduction, and the discrete Fourier transform is the unitary between them.  This is why an FFT grid
  and a plane-wave cutoff are two faces of one fit basis, and why a non-uniform (Becke) grid on a periodic
  cell is still a \f$G=T\f$, delta-family, construction-B basis -- it merely has no \f$\mathbf G\f$-ball dual.

\subsection bst_composition 3.3 Composing the two: the \f$\circ\f$ notation

Most real bases are built by more than one step, written as a **composition** \f$X\circ Y\f$: apply
\f$Y\f$ first, then \f$X\f$ -- read right to left, as in \f$(f\circ g)(x)=f(g(x))\f$.

- Molecular Gaussians are \f$B\circ A\f$: *first* subduce each shell from \f$O(3)\f$ to the site group,
  *then* induce the site up to \f$P\f$.
- GPW is \f$B_T\circ(B_P\circ A)\f$: take the molecular basis and Bloch-sum it over the lattice (a second
  induction, now up to \f$T\f$).  For a full space group \f$T\rtimes P\f$ one induces once more, over the
  little group of \f$\mathbf k\f$.
- Plane waves and \f$Y_{lm}\f$ radial bases are a bare \f$A\f$; deltas and Bloch sums a bare \f$B\f$.

(`doc/BasisSetTaxonomyPlan.md` writes the same pipelines in the same right-to-left convention.)

\section bst_tiers 4. How the code is layered: four tiers, each tagged by one axis

The atom library shows the finished shape and the lattice side now matches it: every IBS flavour is a thin
class templated on an injected **evaluator** (the engine), mixins combine the role aspects, and C++
concepts are the brief spec -- *meet this spec and everything just works*.

| tier | tagged by | what it is | where |
|---|---|---|---|
| **role faces** | role | `Orbital_1E_IBS<T>`, `Orbital_DFT_IBS<T,TFit>`, `FIT_CD_ABS<T>`, ... -- what a Hamiltonian term may ask | `qcBasisSet` (this library) |
| **spec** | G | the concepts an engine must satisfy to drive one G's blocks, plus the evaluator-injected mixins that turn a conforming engine into the role faces | `qchem.BasisSet.Lattice_IBS` (\f$T\f$: `isLattice_{1E,DFT}_Evaluator`, `Lattice::Orbital_1E_IBS<E,T>`); `Radial/IrrepBasisSet.C` (\f$O(3)\f$) |
| **engine** | family | the integral machinery; knows nothing of the spec, merely satisfies it | the `Evaluators/` trees of each engine library |
| **IBS** | (G, family) | one thin class per pair | `PlaneWave_IBS`, `GPW_IBS`, the radial `*_IBS<E>` |

The rules that follow:
- The spec is engine-free, so it lives *below* every engine, in the core.  Its concepts are purely
  structural (no `derived_from`): that is why a Gaussian engine (`GPW_Evaluator`) satisfies the lattice
  spec without inheriting from the plane-wave engine.  The spec is checked where engine meets spec -- a
  `static_assert` beside each concrete IBS -- never inside an engine.
- A concept is named for the G it serves, never for a family.  A family assumption inside a spec is a
  defect.
- A type name reads as the whole taxonomy: `Lattice::Orbital_1E_IBS<GPW_Evaluator, double>` is
  (G = \f$T\f$, role = orbital 1E, family = Gaussian via the GPW engine, scalar = real, i.e. a TRIM block).

\section bst_placement_rules 5. The two placement rules, and the library map

-# **A LIBRARY is an engine.**  A basis library's mass is its integral engine, and an engine factorises over
   the *family*, not over G.  So `UnitCell` inside the Gaussian engine is legitimate -- a lattice sum is a
   Gaussian-engine operation -- and the GPW seam ("GPW is a new evaluator, not a new IBS") stays inside the
   Gaussian library.  A library cut on the G axis is ruled out.
-# **A MODULE carries the G.**  `qchem.BasisSet.Gaussian.Point.*` versus `qchem.BasisSet.Gaussian.Lattice.*`.
   This is a checked invariant: `scripts/audit-basisset-gtags` (ctest `BasisSetGTagAudit`) fails the build
   if a `.Point.` module imports a lattice, or a `.Radial.` module imports a lattice or the molecular
   point-group symmetry.

Corollary -- the **container** tier (the \f$\bigoplus_\Gamma\f$ list, its factory, band structure) is
family-agnostic and becomes its own library only when it must see two or more engines.  True today only for
\f$T\f$ (plane waves and Gaussians; APW will add the radial engine).

| library | module prefix | contents | links |
|---|---|---|---|
| `qcBasisSet` | `qchem.BasisSet.*` | the role faces, the integral caches, `DeltaFit_IBS`, `Projector3`, `Orbital_PP_IBS` (+ its species-field arguments), `GMap`, and the \f$T\f$ spec `Lattice_IBS` | `qcElConfig qcStructure qcMesh` |
| `qcRadial_BS` | `qchem.BasisSet.Radial.*` | the \f$O(3)\f$ engine -- Slater / Gaussian / B-spline radials, the angular integrals, the container + factory.  The one place where the group IS the engine (every integral is radial) | `qcBasisSet` |
| `qcPlaneWave_BS` | `qchem.BasisSet.PlaneWave.*` | the exponential engine: `PW_Evaluator`, the FFT/Poisson grid, `PlaneWave_IBS`, `PlaneWaveFit_IBS`.  A leaf over the core | `qcBasisSet qcSymmetry` |
| `qcGaussian_BS` | `qchem.BasisSet.Gaussian.*` | `Gaussian.Evaluators.*` the engine (no G); `Gaussian.PG_Cart.*` the raw AO block (the seed, see below); `Gaussian.Point.*` \f$G=P\f$: PG_Spherical, PG_LibCint, the SALC container, factory, readers; `Gaussian.Lattice.*` \f$G=T\f$: `LatticeSum1E`, `LatticeScreener`, the spherical lattice view, GPW | `qcBasisSet qcPlaneWave_BS qcSymmetry cint` |
| `qcLattice_BS` | `qchem.BasisSet.Lattice.*` | the \f$T\f$ container: the \f$\bigoplus_{\mathbf k}\f$ `BasisSet` + factories, `BandStructure`, `APW_IBS`, `LAPW_IBS` | `qcPlaneWave_BS qcGaussian_BS qcLASolver` |

**Why `PG_Cart` carries no G tag.**  The Cartesian-Gaussian AO block is the family's raw seed: it answers
`GetAoShells` so that `Point/SymmetryAdaptedBasisSet` can induce it up to \f$P\f$, and it answers
`LatticeSum1E` so that `Lattice/GPW` can Bloch-sum it up to \f$T\f$.  Both constructions act on it from
outside; the block itself carries \f$G=1\f$.  It therefore sits beside `Evaluators/` with no tag, and the
audit does not apply to it.

**Naming.**  `Radial` rather than `Spherical` (which collides with `PG_Spherical`, spherical-harmonic
Gaussians); `Lattice` rather than `Lattice_3D` because \f$T^d\f$ is a row.  `qcSymmetry/{Atom, Molecule,
Lattice_3D}` keep their names: `Symmetry::Molecule` (point groups) versus `Symmetry::Lattice_3D` (space
groups) is a real distinction on the G axis.

\section bst_table 6. The placement table -- everything we have, and everything anticipated

| basis | G | family | construction | role | library |
|---|---|---|---|---|---|
| Slater / Gaussian / B-spline \f$\times\,Y_l\f$, \f$Y_{lm}\f$ | \f$O(3)\f$ | radial | A | Orbital, Fit | `qcRadial_BS` |
| RKB{Slater, Gaussian, B-spline} \f$\times\,\Omega_{\kappa m_j}\f$ | \f$O(3)^*\f$ | radial (RKB) | A | Orbital (DHF) | `qcRadial_BS` |
| PG_Cart / PG_Spherical / PG_LibCint | \f$P\f$ (incl. \f$1\f$) | Gaussian | \f$B\circ A\f$ | Orbital, Fit | `qcGaussian_BS` (`.Point.`) |
| GPW | \f$T\f$ (\f$T\rtimes P\f$ later) | Gaussian | \f$B_T\circ B_P\circ A\f$ | Orbital | `qcGaussian_BS` (`.Lattice.`) |
| plane waves | \f$T\f$ | exponential | A | Orbital | `qcPlaneWave_BS` |
| plane-wave fit | \f$T\f$ | exponential | A | Fit | `qcPlaneWave_BS` |
| uniform delta (\f$v_{xc}\f$ fit) | \f$T\f$ | delta | B (= the PW adjoint) | Fit | `qcBasisSet` (`DeltaFit_IBS`) |
| Becke delta, molecule | \f$1\f$ or \f$P\f$ (site-adapted invariant mesh) | delta | B | Fit | `qcBasisSet` |
| Becke delta, periodic | \f$T\f$ | delta | B | Fit | `qcBasisSet` |
| **anticipated** | | | | | |
| numerical atomic orbitals (SIESTA / FHI-aims) | any | numeric radial \f$\times\,Y_{lm}\f$ | as Gaussians | Orbital | a NEW multi-centre engine library (`qcNumeric_BS`), NOT a `Radial.*` sub-directory |
| APW / LAPW | \f$T\f$ | **composite**: PW outside \f$\otimes\f$ radial\f$\times Y_{lm}\f$ inside | A + sphere matching | Orbital | `qcLattice_BS` (forces the edge to `qcRadial_BS`) |
| molecule in a plane-wave box | \f$T\f$ (large cell) | exponential | A | Orbital | `qcPlaneWave_BS` |
| 2D slabs / 1D polymers | \f$T^2\f$, \f$T^1\f$ (+\f$P\f$) | any | unchanged | -- | unchanged |
| 4-component molecules / SOC crystals | \f$P^*\f$, \f$(T\rtimes P)^*\f$ | RKB Gaussian \f$\times\,\Omega\f$ | \f$B\circ A\f$ | Orbital | `qcGaussian_BS`; the groups are `qcSymmetry` rows |
| Wannier functions | \f$T\f$ (inverse: \f$\mathbf k\to\mathbf R\f$) | derived | \f$B^{-1}\f$ | -- | -- |
| multiwavelets / finite elements | \f$1\f$ or \f$T\f$ | delta-like | B | Fit, Orbital | -- |
| DFT+U / Kleinman-Bylander projectors | -- | radial, subduced | -- | a ROLE (`Orbital_PP_IBS`), not a basis | `qcBasisSet` |

Nothing in the table needs a fourth axis.  The only genuinely new *kind* of entry is the composite family
(APW), and it is what forces the \f$T\f$ container to be its own library.

\section bst_howto 7. Placing a new basis: three questions

-# **Which group labels its blocks?**  That fixes the module tag (`.Point.`, `.Lattice.`, `.Radial.`) and
   the symmetry library it consumes.  If it is a group we do not have (a double group, a space group),
   that is a `qcSymmetry` row first and a basis row second.
-# **Which elementary function does its engine integrate?**  That fixes the LIBRARY.  If the engine exists,
   the new basis is one more thin IBS class in it; if not, it is a new engine library beside the others,
   whatever G it serves.
-# **Is it an orbital basis or a fit basis?**  That fixes the role faces it must implement -- and if it is
   a fit basis, remember that its grid and its family are two independent choices.

The construction axis then follows from the answers to 1 and 2, and the type name of the result reads
back the whole classification.

\section bst_glossary 8. Glossary

- **Branching rule** -- the decomposition of an irrep of a group into irreps of a subgroup
  (\f$l{=}2\to t_{2g}\oplus e_g\f$ under \f$O(3)\downarrow O_h\f$).  The recipe for subduction.
- **Carrier space** -- the vector space on which the matrices of a representation act; for an irrep
  \f$\Gamma\f$, the functions that transform among themselves as \f$\Gamma\f$.  An `IrrepBasisSet` spans
  copies of one.
- **Composition \f$X\circ Y\f$** -- do \f$Y\f$, then \f$X\f$ (right to left).
- **Direct product \f$G_1\times G_2\f$** -- the group of pairs with commuting factors; irreps are pairs of
  irreps.  \f$G_{\rm spatial}\times SU(2)\f$ is the non-relativistic spin structure.
- **Direct sum \f$\oplus\f$** -- concatenation of independent sub-spaces with no matrix elements between
  them; \f$\mathrm{BasisSet}=\bigoplus_\Gamma\mathrm{IrrepBasisSet}_\Gamma\f$.
- **Double group \f$G^*\f$** -- the group obtained by counting a \f$2\pi\f$ rotation as a distinct element,
  needed for half-integer spin; adds the spinor irreps (\f$\Omega_{\kappa m_j}\f$ for the atom).
- **Engine** -- the integral machinery of one family (McMurchie-Davidson, libcint, FFT, radial quadrature).
  A library is an engine.
- **Family** -- a seed type and everything built from it; the axis that fixes the library.
- **G** -- the symmetry group of the Hamiltonian, whose irreps label the blocks.
- **Induction \f$H\uparrow G\f$** (construction B) -- building a representation of \f$G\f$ from a seed local
  to a site with symmetry \f$H\subset G\f$ by carrying it over its orbit and projecting; SALCs, Bloch sums.
- **Irrep** -- an irreducible representation; assumed known.
- **Little group of \f$\mathbf k\f$** -- the point operations that map \f$\mathbf k\f$ to itself (modulo a
  reciprocal-lattice vector); labels space-group irreps together with the star of \f$\mathbf k\f$.
- **\f$O(3)\f$** -- the full rotation-reflection group of three-space about a point; the free atom's \f$G\f$.
- **Orbit** -- the set of images of a site (or function) under all operations of \f$G\f$.
- **Role** -- orbital versus fit; the faces a basis implements.
- **SALC** -- symmetry-adapted linear combination; the result of inducing site-local AOs up to the point group.
- **Seed / analytic seed** -- the elementary function type an engine integrates, before symmetry adaptation.
- **Semidirect product \f$T\rtimes P\f$** -- the space group: translations (normal) with the point group
  acting on them; elements \f$\{R|\mathbf t\}\f$.
- **Site group** -- the subgroup of \f$P\f$ (or of the space group) that leaves a given atom in place.
- **Spec** -- the C++ concepts + mixins that say what an engine must answer to drive one G's blocks.
- **Subduction \f$G_{\rm big}\downarrow G\f$** (construction A) -- restricting a representation of a bigger
  group to a subgroup and re-sorting its carrier into the subgroup's irreps; crystal-field splitting,
  \f$\mathbf k+\mathbf G\f$.
- **\f$T\f$, \f$T^d\f$** -- the lattice translation group in 3 (or \f$d\f$) dimensions; irreps labelled by
  \f$\mathbf k\f$.
- **TRIM** -- a time-reversal-invariant momentum, \f$\mathbf k\equiv-\mathbf k\f$; its Bloch block is real.

Further reading: M. S. Dresselhaus, G. Dresselhaus and A. Jorio, *Group Theory: Application to the Physics
of Condensed Matter* (Springer, 2008) for subduced and induced representations and space groups from a
physicist's side; S. L. Altmann, *Induced Representations in Crystals and Molecules* (Academic Press,
1977) for induction as the unifying construction; C. J. Bradley and A. P. Cracknell, *The Mathematical
Theory of Symmetry in Solids* (Oxford, 1972) for the space-group irreps and double groups.
*/
module;
#include <vector>
#include <memory>

export module qchem.BasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Fit_Types;   // FitQuadrature / VxcFit -- the fit-factory vocabulary
export import qchem.Structure;
export import qchem.Symmetry;
export import qchem.ElectronConfiguration;
export import qchem.Matrix3D;   // Matrix3D -- the crystal point-group ops a periodic basis exposes (IBZ symmetrization)
export import qchem.Symmetry.Lattice_3D.SpaceGroup;   // ReciprocalOp {U|τ} -- the density-symmetrization ops (glide phase)

import qchem.Iterators;
export import qchem.Streamable;

export namespace qchem::BasisSet
{
typedef std::vector<Irrep> irrepv_t;

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
//! \brief A basis set is the direct sum \f$\bigoplus_\Gamma\f$ of one `IrrepBasisSet` per irrep of the
//! Hamiltonian's symmetry group -- see \ref basisset_taxonomy for the road map of every basis type and
//! engine this framework combines, and the definitions behind that sentence.
template <class T> class tBasisSet
    : public virtual Streamable
{
public:
    typedef Orbital_1E_IBS<T> obs_t;

    virtual ~tBasisSet() {};
    virtual size_t   GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;

    virtual FIT_CD_ABS<T>* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams&) const;
    //! \brief The \f$v_{xc}\f$ fit basis for this run: \a fit picks the REPRESENTATION and \a mp the
    //! POINTS, and ONE object comes back carrying both -- because a fit basis is a family of weight vectors
    //! over shared points, so its mesh is constitutive of it, not a sibling return value.
    //!
    //! THIS is the level that chooses between representations (a Bloch BLOCK's own factory builds only the
    //! lineage's fitted basis): \c VxcFit::Delta returns the \f$\delta\f$ basis over the run's XC
    //! quadrature (\c CreateXCQuadrature -- mesh + fold + Shubnikov tags), anything else the lineage's
    //! fitted one.  A molecular caller passes \c Auto and gets its Gaussian auxiliary basis, as always.
    //! \param quad OPTIONAL out: the real-space QUADRATURE the returned basis was built over, when it has
    //! one (\c VxcFit::Delta; default-constructed/empty otherwise).  INJECTION, not a getter (user,
    //! 2026-08-23): this factory CREATES the quadrature, so handing it to a second collaborator is what a
    //! creator does -- whereas an accessor on the fit basis would hand out its own private state, which is
    //! exactly the escape the 2026-08-22 arc spent itself closing.
    //!
    //! The WHOLE bundle since 2026-08-24, not just its mesh.  Two of its fields were reaching consumers by
    //! the other route -- as \c FIT_SF_ABS::Symmetrize / \c SymmetrizeSpin, i.e. as OPERATIONS on the fit
    //! face, which is precisely where the orbit fold and the Shubnikov tags do not belong (a fit basis is a
    //! family of functions; star-averaging a coefficient vector over a crystal group is not a fitting
    //! question).  Handing the sibling fields the same way as the mesh removes both members and lets the
    //! consumer call the free \c SymmetrizeValues / \c SymmetrizeValuesSigned directly.  The mesh alone is
    //! still what the ATOMIC PARTITION observable needs, which is why it was the first field to travel.
    virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams&,
                                                VxcFit fit=VxcFit::Auto,
                                                FitQuadrature* quad=nullptr) const;
    //! The DELTA-fit sibling of \c CreateVxcFitBasisSet: the finished real-space XC quadrature (mesh +
    //! symmetry fold), assembled by the basis -- which owns the cell and the §3-imposed ops -- so the
    //! Hamiltonian does no mesh work.  Default/plain path: the Structure's integration mesh, no fold.
    virtual FitQuadrature   CreateXCQuadrature (const Structure* cl, const qcMesh::MeshParams&) const;

    //! The crystal RECIPROCAL point group as \f${U|\tau}\f$ ops for IBZ density symmetrization.  Default {} =
    //! trivial {E} = no-op -- molecules / Γ / unfolded bases.  A periodic GPW basis returns the ops WHEN it
    //! folds the mesh (imposeSymmetry), so the composite density ctor-injects them and the reduced density is star-
    //! averaged with the glide phase \f$e^{+2\pi i(Um)\cdot\tau}\f$ (doc/GPWPlan1.md items 3 + 5).  It is a basis
    //! property (the basis computes the space group), so the density/WF read it here rather than via a setter.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetReciprocalPointOps() const {return {};}

    //! The FULL DETECTED crystal reciprocal point group \f${U|\tau}\f$ -- populated on every periodic run,
    //! IMPOSED or not (unlike \c GetReciprocalPointOps, which is the ops-to-impose set and stays empty on a
    //! free run).  This is what the §3 order-parameter diagnostic (doc/SymmetryUpgradePlan.md) measures the
    //! converged density against: a free run reports the symmetry it actually found, per op
    //! (\c SymmetryDefects).  Default {} = no lattice / nothing detected.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetDetectedReciprocalOps() const {return {};}

    // Iterate() with no type argument yields the base obs_t* directly (no cast);
    // Iterate<D>() dynamic_cast's each IBS to the requested derived type D.
    // Built on the two primitives below, so storage stays private to the
    // concrete BasisSet.
    auto Iterate() const {return IndexProxy<const tBasisSet>(this, GetNumIBS());}
    template <class D> auto Iterate() const
    {
        return D_IndexProxy<const D, const tBasisSet>(this, GetNumIBS());
    }
    const obs_t* operator[](size_t i) const {return GetIBS(i);}

    virtual size_t GetNumIBS() const=0; //!< Number of Irrep basis sets (public since Step 3c-2: the
                                        //!< mixed-aware per-index walk needs the count beside GetRealIBS).
    //! \brief The CROSS-SCALAR block view (doc/RealComplexPlan.md Step 3c-2): non-null iff block \a i is a
    //! REAL block inside a complex-faced set (a TRIM block after the Step-3c-3 factory decision).  The
    //! same-scalar face stays \c GetIBS / \c Iterate; a mixed-aware consumer (MakeIrrepWFs, the GPW
    //! preflight) checks this FIRST and falls back to \c GetIBS.  Default null: a homogeneous set -- and
    //! every real-faced set, whose blocks are already real via \c GetIBS -- has no cross-scalar children.
    virtual const Orbital_1E_IBS<double>* GetRealIBS(size_t) const {return nullptr;}

protected:
    virtual const obs_t* GetIBS(size_t) const=0; //The only storage-specific primitive.
};

typedef tBasisSet<double>    Real_BS;
typedef tBasisSet<dcmplx> Complex_BS;

}//namespace