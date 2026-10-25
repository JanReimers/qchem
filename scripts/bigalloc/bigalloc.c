// bigalloc: LD_PRELOAD shim -- log every allocation >= BIGALLOC_MIN bytes (default 1 MB) with a backtrace, and every
// free of a logged pointer, to $BIGALLOC_OUT.  `live.py` turns the log into the LIVE SET AT PEAK, grouped by qchem
// source line (addr2line, inlines expanded).  Built by `scripts/bigalloc/run` on first use.  (No heaptrack/valgrind
// on this box; written 2026-10-25 for the plane-wave footprint work, doc/OpenWork.md "Dense PW cost".)
#define _GNU_SOURCE
#include <dlfcn.h>
#include <execinfo.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#define write(a,b,c) (void)!write(a,b,c)
static void* (*r_malloc)(size_t);
static void  (*r_free)(void*);
static void* (*r_calloc)(size_t,size_t);
static void* (*r_realloc)(void*,size_t);
static int   (*r_posix_memalign)(void**,size_t,size_t);
static void* (*r_aligned_alloc)(size_t,size_t);
static int fd=-1; static size_t minsz=1<<20;
static __thread int busy=0;
static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;
#define NSLOT (1<<18)
static void* slots[NSLOT];   // logged pointers (open addressing; (void*)1 = tombstone)

static size_t h(void* p){ return ((size_t)p>>4)*2654435761u & (NSLOT-1); }
static void put(void* p){ size_t i=h(p); while(slots[i] && slots[i]!=(void*)1) i=(i+1)&(NSLOT-1); slots[i]=p; }
static int  take(void* p){ size_t i=h(p); while(slots[i]){ if(slots[i]==p){ slots[i]=(void*)1; return 1;} i=(i+1)&(NSLOT-1);} return 0; }

// dlsym itself allocates: serve those early calls from a static arena.
static char boot[1<<16]; static size_t bootn=0; static int inited=0, initing=0;
static void resolve(void)
{
    if (inited||initing) return; initing=1;
    r_malloc=dlsym(RTLD_NEXT,"malloc"); r_free=dlsym(RTLD_NEXT,"free"); r_calloc=dlsym(RTLD_NEXT,"calloc");
    r_realloc=dlsym(RTLD_NEXT,"realloc"); r_posix_memalign=dlsym(RTLD_NEXT,"posix_memalign"); r_aligned_alloc=dlsym(RTLD_NEXT,"aligned_alloc");
    inited=1; initing=0;
}
__attribute__((constructor)) static void init(void)
{
    resolve();
    const char* o=getenv("BIGALLOC_OUT"); const char* m=getenv("BIGALLOC_MIN");
    if (m) minsz=strtoull(m,0,10);
    if (o) fd=open(o,O_WRONLY|O_CREAT|O_TRUNC,0644);
    void* b[2]; busy=1; backtrace(b,2); busy=0;   // warm up the libgcc unwinder outside malloc
}
static void logalloc(void* p, size_t n)
{
    if (fd<0||!p||n<minsz||busy) return;
    busy=1; pthread_mutex_lock(&mu);
    char line[64]; int k=snprintf(line,sizeof line,"A %p %zu\n",p,n); write(fd,line,k);
    void* bt[40]; int d=backtrace(bt,40); backtrace_symbols_fd(bt+2,d-2,fd); write(fd,"E\n",2);
    put(p); pthread_mutex_unlock(&mu); busy=0;
}
static void logfree(void* p)
{
    if (fd<0||!p||busy) return;
    pthread_mutex_lock(&mu);
    if (take(p)) { char line[40]; int k=snprintf(line,sizeof line,"F %p\n",p); write(fd,line,k); }
    pthread_mutex_unlock(&mu);
}
void* malloc(size_t n){ resolve(); if(!r_malloc){ void* p=boot+bootn; bootn+=(n+15)&~15ul; return p; } void* p=r_malloc(n); logalloc(p,n); return p; }
void  free(void* p){ if((char*)p>=boot && (char*)p<boot+sizeof boot) return; logfree(p); if(r_free) r_free(p); }
void* calloc(size_t a,size_t b){ resolve(); if(!r_calloc){ void* p=boot+bootn; bootn+=(a*b+15)&~15ul; memset(p,0,a*b); return p; } void* p=r_calloc(a,b); logalloc(p,a*b); return p; }
void* realloc(void* q,size_t n){ resolve(); logfree(q); void* p=r_realloc(q,n); logalloc(p,n); return p; }
int posix_memalign(void** pp,size_t a,size_t n){ resolve(); int r=r_posix_memalign(pp,a,n); if(!r) logalloc(*pp,n); return r; }
void* aligned_alloc(size_t a,size_t n){ resolve(); void* p=r_aligned_alloc(a,n); logalloc(p,n); return p; }
