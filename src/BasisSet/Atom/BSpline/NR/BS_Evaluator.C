// File: BasisSet/Atom/BSpline/NR/BS_Evaluator.C
module;
export module BasisSet.Atom.BSpline.NR.BS_Evaluator;
export import qchem.BasisSet.Atom.Gaussian.NR.BS_Evaluator;
import qchem.BasisSet.Atom.BSpline.SplineGrouper;
import qchem.BasisSet.Atom.BSpline.Rk;
import qchem.Basisset.Atom.BSpline.GLQuadrature;

export template <size_t K> class BSpline_BS 
    : public virtual BS_Evaluator
{
public:
    virtual void Register(IBS_Evaluator *);
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual RVec loop_4_direct  (size_t id, size_t la, size_t lc) const;
    virtual RVec loop_4_exchange(size_t id, size_t la, size_t lc) const;
protected:
    void BuildCache(size_t lmax);
private:
    const GLCache* GetGL(size_t l) const;
    SplineGrouper<K> grouper;
    BSpline::RkCache<K>* itsRkCache;
};

