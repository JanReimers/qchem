#include <cstdio>
struct Dyn  { virtual ~Dyn()=default; };                    // tDynamic_HT: the abstract face (the cast's source)
struct Res  { virtual ~Res()=default; };                    // tRealSpaceResidual: the target face, unrelated
struct Imp  : public virtual Dyn {};                        // tDynamic_HT_Imp: shares the virtual base
struct Bad  : public virtual Dyn, public virtual Res, protected Imp {};   // face BEFORE the protected base
struct Good : public virtual Dyn, protected Imp, public virtual Res {};   // face AFTER it
template<class T> void t(const char* n){ T x; Dyn* d=&x; Res* r=dynamic_cast<Res*>(d); std::printf("%-5s Dyn->Res %s\n",n,r?"ok":"NULL"); }
int main(){ t<Bad>("Bad"); t<Good>("Good"); }
