// cmake/OpenMPRuntimeCheck.cpp  Configure's proof that a program linked the tree's way LOADS the libomp it was linked against.
// Built and run by the root CMakeLists.txt with the tree's OpenMP flags, OpenMP_omp_LIBRARY and CMAKE_BUILD_RPATH.  It prints the
// canonical path of the shared object that holds omp_get_max_threads at RUN time; configure compares it with the linked library.
// Earned 2026-11-23 (lane G, clang 23): the tree linked /opt's libomp, RUNPATH lacked its directory, and every binary silently
// loaded apt's older /usr/lib libomp instead.
#include <dlfcn.h>
#include <climits>
#include <cstdio>
#include <cstdlib>
extern "C" int omp_get_max_threads();
int main()
{
    Dl_info info{};
    if (omp_get_max_threads()<1 || !dladdr(reinterpret_cast<const void*>(&omp_get_max_threads), &info) || !info.dli_fname) return 1;
    char real[PATH_MAX];
    std::printf("%s", realpath(info.dli_fname, real) ? real : info.dli_fname);
    return 0;
}
