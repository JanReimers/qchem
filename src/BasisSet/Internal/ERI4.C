// File ERI4.C  Symmetric containters for and ERI (Electron Repulsion Integral) 4 index super matrix.
module;
export module qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI4T;
export import qchem.Types;

export class ERI4 : public ERI4T<double,smat_t>
{
public:
    typedef ERI4T<double,smat_t> Base;
    ERI4() : Base() {};
    ERI4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    friend void MatMul(rsmat_t& Sab, const ERI4& gabcd,const rsmat_t& Scd);
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

export bool operator==(const ERI4& a, const ERI4& b);


