// File: MolecularMesh.C  Becke-partitioned multi-center molecular integration mesh.
//
// Geometry-AWARE consumer of qcMesh (it needs atom positions), so it lives in qcStructure, not in
// qcMesh.  Clean-room replacement for src/Structure/Internal/MoleculeMesh.C: a from-scratch Becke
// rewrite that, crucially, handles the coincident-atom degeneracy the old code divided-by-zero on.
//
//   for each atom a:  Mesh m = ProductMesh(radial, angular);  m.ShiftOrigin(R_a);
//                     reweight each point by the Becke cell function  p_a(r) / sum_b p_b(r)
//                     append to the combined mesh
//
// (Atom-size adjustment of the cell boundaries is NOT applied -- matching the old mesh, so the
// pinned molecular-DFT energies do not move when the consumers migrate.)
module;
export module qchem.Structure.MolecularMesh;
export import qchem.Mesh;          // qcMesh::Mesh, qcMesh::MeshParams (the result + the knobs)
import qchem.Structure;             // Structure (brings the OLD global Mesh -- we always qualify qcMesh::)

//! \brief Becke fuzzy-Voronoi integration mesh for a finite structure.  mp.beckeOrder is the number
//! of smoothing iterations of Becke's cell polynomial (Becke 1988 recommends 3).  Coincident atoms
//! (R_ab = 0) are handled explicitly: their pair coordinate mu is taken as 0, so each contributes a
//! half-weight cell and a coincident dimer integrates to exactly the single-atom result.
export qcMesh::Mesh MakeMolecularMesh(const Structure& cl, const qcMesh::MeshParams& mp);
