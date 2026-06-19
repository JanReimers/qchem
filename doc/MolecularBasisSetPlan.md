# Molecular basis set upgrade

## High level goals, big picture (we may have to re-order these): 

    
1. Use cache4 in src/BasisSet/Internal/Cache4.C, for Molecular integrals.  This is a way to store all information
    specific to a four index combination of basis functions.  Just derive from Cacheable4 and put any data memebers
    you want in the derived class.  Used in 4 index loops for HF Direct/Exchange integrals.
2. Use the Evalautor idea for PG basis sets, move src/BasisSet/Atom/Evaluators/Evaluator.C (without Getl() and the Angular class) out to src/BasisSet/Internal/Evaluator.C.
3. Move 1E Integrals_Overlap,Integrals_Kinetic,Integrals_Nuclear (src/BasisSet/Atom/IrrepBasisSet.C) out to src/BasisSet/Orbital_1E_IBS.C, and so on.  Basiscally all those i,j loops to build matrices should be basis set agnostic (even for solids yes?).
4. Move ERI loops in Orbital_HF_IBS<E>::MakeDirect/MakeExchange out to src/BasisSet/Orbital_HF_IBS.C.  This might be challenging becuase of Ak calculations for atoms.  We may need a new abstraction here.
5. Once we have a very clean PG implementation with as much generic code as possible at then BasisSet level (out of Molecule and Atom subfolders),
we should add more molecular basis set types.  If we did everything right this just means definining new Evaluators: 

    - SphericalGaussian using M&D
    - PolarizedGausian using https://github.com/sunqm/libcint and/or (https://github.com/ValeevGroup/libintx/) instead of M&D
    - SphericalGaussian using libCint and/or libintX M&D
    - More modern algorithms.  For example ACE method TAKESHI YANAI et. al. International Journal of Quantum Chemistry, Vol. 76, 396᎐406 (2000).

6. Define a BasisSetSourse interface with implementations for
    - read basis sets (Start with GAUSSIAN94 see src/BasisSet/Molecule/PolarizedGaussian/Reader.C) from the local file system.
    - Pull basis sets from https://www.basissetexchange.org/ and/or https://molssi-bse.github.io/basis_set_exchange/
    
- But before do any of that I think the current PG implementation needs a major refactor as described in the next section.

## Stage 1: PG refactor 
- Currently the PG code uses 4 way virtual dispatch mechanism to handle combinations of contracted radial functions and primative gaussians in integrals.  Unlike julia, c++ has no native support for multiple dispatch so the code is complex and possibly slow.  Lots of dynamic casts.  I now realize this was too clever by half, over use of the dependency inverson  principle.
- I think a better approach is to just have one type of radial function with contraction coefficients.  Primatives just have one coefficient set to 1.0 or the normalization contstant.
- If I was doing this myself I would create a now tree src/BasisSet/Molecule/PolarizedGaussian1 (namespace BasisSet::Molecule::PolarizedGaussian1), that temporarily uses src/BasisSet/Molecule/PolarizedGaussian/MnD and src/BasisSet/Molecule/PolarizedGaussian/Readers, polarizations, Reader.C, Symmetry.C  etc as is.
- And build up the new (simpler) radial function mechanism from scratch.
- But I always want to hear your opinion on the best way to proceed.
- Once it is proven working for H20, N2 HF we scrub the old PG, rename the new one and move on the stage2.


## Stage 2: New PG Use cache_4 
- Following McMurchie & Davidson (M&D) the code uses the concept of a charge distribution Ω_ab (GaussianCD in code src/BasisSet/Molecule/PolarizedGaussian/Internal/CDCache.C) for storing any data associated and ab (primative) basis function pair.
- Right now these are cached inside the PG classes.  This cacheing should move to the src/BasisSet/Internal/DB_Cache.C system.  
- in src/BasisSet/Internal we need to make a Cache2.C which supports 2 index caching and maybe looping.
- Then the new class PG_cache4 : public virtual Cacheable {...};  class can store pointers to Ω_ab an Ω_cd and any thing else it needs.
- We can see ehat else is needed in PG_cache4 by looking in src/BasisSet/Molecule/PolarizedGaussian/Radial/Imp/GaussianRF.C GaussianRF::Integrate4C.   I just see double lambda=2*Pi52/(ab.AlphaP*cd.AlphaP*sqrt(ab.AlphaP+cd.AlphaP)); //M&D 3.31 as a candidate.

## Stage 3: New PG Use cache_3 for DFT
- I think this is self evident.
- Once this is done we can pause and think.  



