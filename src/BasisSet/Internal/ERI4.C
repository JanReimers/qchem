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
    friend rsmat_t MatMul(const ERI4& gabcd,const rsmat_t& Scd);
    friend rsmat_t MatMul(const rsmat_t& Sab, const ERI4& gabcd);
private:
    static double contract(const rsmat_t& A,const rsmat_t& B);
};

export class M4 : public ERI4T<double,mat_t>
{
public:
    typedef ERI4T<double,mat_t> Base;
    M4() : Base() {};
    M4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    friend rmat_t MatMul(const M4& gabcd,const rmat_t& Scd);
    friend rmat_t MatMul(const rmat_t& Sab, const ERI4& gabcd);
private:
    static double contract(const rmat_t& A,const rmat_t& B);
};
 


