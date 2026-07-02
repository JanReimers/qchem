// File ERI4.C  Symmetric containters for and ERI (Electron Repulsion Integral) 4 index super matrix.
module;
export module qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI4T;
export import qchem.Types;

namespace qchem {

export class ERI4 : public ERI4T<double,smat_t>
{
public:
    typedef ERI4T<double,smat_t> Base;
    ERI4() : Base() {};
    ERI4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    //! Contract this block against a cd-density: Sab += sum_cd J(a,b) ⊙ Dcd (each inner block summed to a
    //! scalar -- the fast, localized direction).  Use ScatterBoth() to also feed the bra-ket partner.
    void MatMul(rsmat_t& Sab, const rsmat_t& Dcd) const;
    //! Fused bra-ket scatter of this canonical (i,j) block into BOTH Fock sub-blocks in ONE pass over J:
    //!   Si += J·Dj    (inner block ⊙ Dj -> scalar; localized)
    //!   Sj += Jᵀ·Di   (scalar Di(a,b) times the WHOLE inner block, added to Sj -- never a transposed
    //!                  gather, so no 4× scatter penalty)
    //! Equivalent to MatMul(Si,Dj) plus Transpose().MatMul(Sj,Di) but reads J once and needs only the one
    //! canonical block in the cache.  All the (a,b) symmetry bookkeeping is here.  See doc/ERI4Rework.md.
    void ScatterBoth(rsmat_t& Si, rsmat_t& Sj, const rsmat_t& Di, const rsmat_t& Dj) const;
    friend void MatMul(rsmat_t& Sab, const ERI4& gabcd,const rsmat_t& Scd); //thin delegate to the member
    ERI4 Transpose() const;  //convert Jabcd->Jcdab;
};

export class M4 : public ERI4T<double,mat_t>
{
public:
    typedef ERI4T<double,mat_t> Base;
    M4() : Base() {};
    M4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    friend rmat_t MatMul(const M4& gabcd,const rmat_t& Scd);
};

//
//  The may be usefull for unit testing.
// 
export bool operator==(const ERI4& a, const ERI4& b);
export double fnorm(const ERI4& a, const ERI4& b);
export double relative_fnorm(const ERI4& a, const ERI4& b);


} // namespace qchem