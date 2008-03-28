/**************************************************************************
Isomalloc:
  A way to allocate memory at the same address on every processor.
This enables linked data structures, like thread stacks, to be migrated
to the same address range on other processors.  This is similar to an
explicitly managed shared memory system.

  The memory is used and released via the mmap()/mumap() calls, so unused
memory does not take any (RAM, swap or disk) space.

  The way it's implemented is that each processor claims some section 
of the available virtual address space, and satisfies all new allocations
from that space.  Migrating structures use whatever space they started with.

Written for migratable threads by Milind Bhandarkar around August 2000;
generalized by Orion Lawlor November 2001.
*/
#include "converse.h"
#include "memory-isomalloc.h"

#define CMK_THREADS_DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <errno.h> /* just so I can find dynamically-linked symbols */

static int read_randomflag(void)
{
  FILE *fp;
  int random_flag;
  CmiLock(smp_mutex);
  fp = fopen("/proc/sys/kernel/randomize_va_space", "r");
  if (fp != NULL) {
    fscanf(fp, "%d", &random_flag);
  }
  else {
    random_flag = -1;
  }
  CmiUnlock(smp_mutex);
  return random_flag;
}

struct CmiIsomallocBlock {
      CmiInt8 slot; /*First mapped slot*/
      CmiInt8 length; /*Length of (user portion of) mapping, in bytes*/
};
typedef struct CmiIsomallocBlock CmiIsomallocBlock;

/*Convert a heap block pointer to/from a CmiIsomallocBlock header*/
static void *block2pointer(CmiIsomallocBlock *blockHeader) {
	return (void *)(blockHeader+1);
}
static CmiIsomallocBlock *pointer2block(void *heapBlock) {
	return ((CmiIsomallocBlock *)heapBlock)-1;
}

/*Integral type to be used for pointer arithmetic:*/
typedef size_t memRange_t;

/*Size in bytes of a single slot*/
static size_t slotsize;

/*Total number of slots per processor*/
static CmiInt8 numslots=0;

/*Start and end of isomalloc-managed addresses.
If isomallocStart==NULL, isomalloc is disabled.
*/
static char *isomallocStart=NULL;
static char *isomallocEnd=NULL;

/*Utility conversion functions*/
static void *slot2addr(CmiInt8 slot) {
	return isomallocStart+((memRange_t)slotsize)*((memRange_t)slot);
}
static int slot2pe(CmiInt8 slot) {
	return slot/numslots;
}
static CmiInt8 pe2slot(int pe) {
	return pe*numslots;
}
/*Return the number of slots in a block with n user data bytes*/
static int length2slots(int nBytes) {
	return (sizeof(CmiIsomallocBlock)+nBytes+slotsize-1)/slotsize;
}

typedef struct _slotblock
{
  CmiInt8 startslot;
  CmiInt8 nslots;
} slotblock;

typedef struct _slotset
{
  int maxbuf;
  slotblock *buf;
  CmiInt8 emptyslots;
} slotset;

/*
 * creates a new slotset of nslots entries, starting with all
 * empty slots. The slot numbers are [startslot,startslot+nslot-1]
 */
static slotset *
new_slotset(CmiInt8 startslot, CmiInt8 nslots)
{
  int i;
  slotset *ss = (slotset*) malloc_reentrant(sizeof(slotset));
  _MEMCHECK(ss);
  ss->maxbuf = 16;
  ss->buf = (slotblock *) malloc_reentrant(sizeof(slotblock)*ss->maxbuf);
  _MEMCHECK(ss->buf);
  ss->emptyslots = nslots;
  ss->buf[0].startslot = startslot;
  ss->buf[0].nslots = nslots;
  for (i=1; i<ss->maxbuf; i++)
    ss->buf[i].nslots = 0;
  return ss;
}

/*
 * returns new block of empty slots. if it cannot find any, returns (-1).
 */
static CmiInt8
get_slots(slotset *ss, CmiInt8 nslots)
{
  int i;
  if(ss->emptyslots < nslots)
    return (-1);
  for(i=0;i<(ss->maxbuf);i++)
    if(ss->buf[i].nslots >= nslots)
      return ss->buf[i].startslot;
  return (-1);
}

/* just adds a slotblock to an empty position in the given slotset. */
static void
add_slots(slotset *ss, CmiInt8 sslot, CmiInt8 nslots)
{
  int pos; 
  int emptypos = -1;
  if (nslots == 0)
    return;
  for (pos=0; pos < (ss->maxbuf); pos++) {
    if (ss->buf[pos].nslots == 0) {
      emptypos = pos;
      break; /* found empty slotblock */
    }
  }
  if (emptypos == (-1)) /*no empty slotblock found */
  {
    int i;
    int newsize = ss->maxbuf*2;
    slotblock *newbuf = (slotblock *) malloc_reentrant(sizeof(slotblock)*newsize);
    _MEMCHECK(newbuf);
    for (i=0; i<(ss->maxbuf); i++)
      newbuf[i] = ss->buf[i];
    for (i=ss->maxbuf; i<newsize; i++)
      newbuf[i].nslots  = 0;
    free_reentrant(ss->buf);
    ss->buf = newbuf;
    emptypos = ss->maxbuf;
    ss->maxbuf = newsize;
  }
  ss->buf[emptypos].startslot = sslot;
  ss->buf[emptypos].nslots = nslots;
  ss->emptyslots += nslots;
  return;
}

/* grab a slotblock with specified range of blocks
 * this is different from get_slots, since it pre-specifies the
 * slots to be grabbed.
 */
static void
grab_slots(slotset *ss, CmiInt8 sslot, CmiInt8 nslots)
{
  CmiInt8 pos, eslot, e;
  eslot = sslot + nslots;
  for (pos=0; pos < (ss->maxbuf); pos++)
  {
    if (ss->buf[pos].nslots == 0)
      continue;
    e = ss->buf[pos].startslot + ss->buf[pos].nslots;
    if(sslot >= ss->buf[pos].startslot && eslot <= e)
    {
      CmiInt8 old_nslots;
      old_nslots = ss->buf[pos].nslots;
      ss->buf[pos].nslots = sslot - ss->buf[pos].startslot;
      ss->emptyslots -= (old_nslots - ss->buf[pos].nslots);
      add_slots(ss, sslot + nslots, old_nslots - ss->buf[pos].nslots - nslots);
      return;
    }
  }
  CmiAbort("requested a non-existent slotblock\n");
}

/*
 * Frees slot by adding it to one of the blocks of empty slots.
 * this slotblock is one which is contiguous with the slots to be freed.
 * if it cannot find such a slotblock, it creates a new slotblock.
 * If the buffer fills up, it adds up extra buffer space.
 */
static void
free_slots(slotset *ss, CmiInt8 sslot, CmiInt8 nslots)
{
  int pos;
  /* eslot is the ending slot of the block to be freed */
  CmiInt8 eslot = sslot + nslots;
  for (pos=0; pos < (ss->maxbuf); pos++)
  {
    CmiInt8 e = ss->buf[pos].startslot + ss->buf[pos].nslots;
    if (ss->buf[pos].nslots == 0)
      continue;
    /* e is the ending slot of pos'th slotblock */
    if (e == sslot) /* append to the current slotblock */
    {
	    ss->buf[pos].nslots += nslots;
	    ss->emptyslots += nslots;
	    return;
    }
    if(eslot == ss->buf[pos].startslot) /* prepend to the current slotblock */
    {
	    ss->buf[pos].startslot = sslot;
	    ss->buf[pos].nslots += nslots;
	    ss->emptyslots += nslots;
	    return;
    }
  }
  /* if we are here, it means we could not find a slotblock that the */
  /* block to be freed was combined with. */
  add_slots(ss, sslot, nslots);
}

/*
 * destroys slotset-- currently unused
static void
delete_slotset(slotset* ss)
{
  free_reentrant(ss->buf);
  free_reentrant(ss);
}
 */

#if CMK_THREADS_DEBUG
static void
print_slots(slotset *ss)
{
  int i;
  CmiPrintf("[%d] maxbuf = %d\n", CmiMyPe(), ss->maxbuf);
  CmiPrintf("[%d] emptyslots = %d\n", CmiMyPe(), ss->emptyslots);
  for(i=0;i<ss->maxbuf;i++) {
    if(ss->buf[i].nslots)
      CmiPrintf("[%d] (%d, %d) \n", CmiMyPe(), ss->buf[i].startslot, 
          ss->buf[i].nslots);
  }
}
#else
#  define print_slots(ss) /*empty*/
#endif

/*This version of the allocate/deallocate calls are used if the 
real mmap versions are disabled.*/
static int disabled_map_warned=0;
static void *disabled_map(int nBytes) 
{
	if (!disabled_map_warned) {
		disabled_map_warned=1;
		if (CmiMyPe()==0)
			CmiError("charm isomalloc.c> Warning: since mmap() doesn't work,"
			" you won't be able to migrate threads\n");
	}
	return malloc(nBytes);
}
static void disabled_unmap(void *bk) {
	free(bk);
}

/*Turn off isomalloc memory, for the given reason*/
static void disable_isomalloc(const char *why)
{
    isomallocStart=NULL;
    isomallocEnd=NULL;
#if CMK_THREADS_DEBUG
    CmiPrintf("[%d] isomalloc.c> Disabling isomalloc because %s\n",CmiMyPe(),why);
#endif
}

#if ! CMK_HAS_MMAP
/****************** Manipulate memory map (Win32 non-version) *****************/
static void *call_mmap_fixed(void *addr,size_t len) {
	CmiAbort("isomalloc.c: mmap_fixed should never be called here.");
}
static void *call_mmap_anywhere(size_t len) {
	CmiAbort("isomalloc.c: mmap_anywhere should never be called here.");
}
static void call_munmap(void *addr,size_t len) {
	CmiAbort("isomalloc.c: munmap should never be called here.");
}

static int 
init_map(char **argv)
{
  return 0; /*Isomalloc never works without mmap*/
}
#else /* CMK_HAS_MMAP */
/****************** Manipulate memory map (UNIX version) *****************/
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#if !CMK_HAS_MMAP_ANON
CpvStaticDeclare(int, zerofd); /*File descriptor for /dev/zero, for mmap*/
#endif

/**
 * Maps this address with these flags.
 */
static void *call_mmap(void *addr,size_t len, int flags) {
  void *ret=mmap(addr, len, PROT_READ|PROT_WRITE,
#if CMK_HAS_MMAP_ANON
	    flags|MAP_PRIVATE|MAP_ANON,-1,
#else
	    flags|MAP_PRIVATE,CpvAccess(zerofd),
#endif
	    0);
  if (ret==((void*)(-1))) return (void *)0; /* all-ones means failure */
  else return ret;
}
static void *call_mmap_fixed(void *addr,size_t len) {
	return call_mmap(addr,len,MAP_FIXED);
}
static void *call_mmap_anywhere(size_t len) {
	return call_mmap((void *)0,len,0);
}

/* Unmaps this address range */
static void call_munmap(void *addr,size_t len) {
  int retval;
  if (addr == 0) return; /* NULL address is never mapped */ 
  retval = munmap(addr, len);
  if (retval==(-1))
    CmiAbort("munmap call failed to deallocate requested memory.\n");
}

static int 
init_map(char **argv)
{
#if CMK_HAS_MMAP_ANON
  /*Don't need /dev/zero*/
#else
  CpvInitialize(int, zerofd);  
  CpvAccess(zerofd) = open("/dev/zero", O_RDWR);
  if(CpvAccess(zerofd)<0)
    return 0; /* Cannot open /dev/zero or use MMAP_ANON, so can't mmap memory */
#endif
  return 1;
}

#endif /* UNIX memory map */


/**
 * maps the virtual memory associated with slot using mmap
 */
static CmiIsomallocBlock *
map_slots(CmiInt8 slot, CmiInt8 nslots)
{
  void *pa;
  void *addr=slot2addr(slot);
  pa = call_mmap_fixed(addr, slotsize*nslots);
  
  if (pa==NULL) 
  { /*Map just failed completely*/
#if CMK_THREADS_DEBUG
    perror("mmap failed");
    CmiPrintf("[%d] tried to mmap %p, but encountered error\n",CmiMyPe(),addr);
#endif
    return NULL;
  }
  if (pa != addr)
  { /*Map worked, but gave us back the wrong place*/
#if CMK_THREADS_DEBUG
    CmiPrintf("[%d] tried to mmap %p, but got %p back\n",CmiMyPe(),addr,pa);
#endif
    call_munmap(addr,slotsize*nslots);
    return NULL;
  }
#if CMK_THREADS_DEBUG
  CmiPrintf("[%d] mmap'd slots %d-%d to address %p\n",CmiMyPe(),
	    slot,slot+nslots-1,addr);
#endif
  return (CmiIsomallocBlock *)pa;
}

/*
 * unmaps the virtual memory associated with slot using munmap
 */
static void
unmap_slots(CmiInt8 slot, CmiInt8 nslots)
{
  void *addr=slot2addr(slot);
  call_munmap(addr, slotsize*nslots);
#if CMK_THREADS_DEBUG
  CmiPrintf("[%d] munmap'd slots %d-%d from address %p\n",CmiMyPe(),
	    slot,slot+nslots-1,addr);
#endif
}

static void map_failed(CmiInt8 s,CmiInt8 n)
{
  void *addr=slot2addr(s);
  CmiError("charm isomalloc.c> map failed to allocate %d bytes at %p.\n",
      slotsize*n, addr);
  CmiAbort("Exiting\n");
}



/************ Address space voodoo: find free address range **********/

CpvStaticDeclare(slotset *, myss); /*My managed slots*/

/*This struct describes a range of virtual addresses*/
typedef struct {
  char *start; /*First byte of region*/
  memRange_t len; /*Number of bytes in region*/
  const char *type; /*String describing memory in region (debugging only)*/
} memRegion_t;

/*Estimate the top of the current stack*/
static void *__cur_stack_frame(void)
{
  char __dummy;
  void *top_of_stack=(void *)&__dummy;
  return top_of_stack;
}
/*Estimate the location of the static data region*/
static void *__static_data_loc(void)
{
  static char __dummy;
  return (void *)&__dummy;
}

/*Pointer comparison is in these subroutines, because
  comparing arbitrary pointers is nonportable and tricky.
*/
static int pointer_lt(const char *a,const char *b) {
	return ((memRange_t)a)<((memRange_t)b);
}
static int pointer_ge(const char *a,const char *b) {
	return ((memRange_t)a)>=((memRange_t)b);
}

static char *pmin(char *a,char *b) {return pointer_lt(a,b)?a:b;}
static char *pmax(char *a,char *b) {return pointer_lt(a,b)?b:a;}

const static memRange_t meg=1024u*1024u; /*One megabyte*/
const static memRange_t gig=1024u*1024u*1024u; /*One gigabyte*/

/*Check if this memory location is usable.  
  If not, return 1.
*/
static int bad_location(char *loc) {
  void *addr;
  addr=call_mmap_fixed(loc,slotsize);
  if (addr==NULL) {
#if CMK_THREADS_DEBUG
    CmiPrintf("[%d] Skipping unmappable space at %p\n",CmiMyPe(),loc);
#endif
    return 1; /*No good*/
  }
  call_munmap(addr,slotsize);
  return 0; /*This works*/
}

/* Split this range up into n pieces, returning the size of each piece */
static memRange_t divide_range(memRange_t len,int n) {
	return (len+1)/n;
}

/* Return if this memory region has *any* good parts. */
static int partially_good(char *start,memRange_t len,int n) {
  int i;
  memRange_t quant=divide_range(len,n);
  for (i=0;i<n;i++)
    if (!bad_location(start+i*quant)) return 1; /* it's got some good parts */
  return 0; /* all locations are bad */
}

/* Return if this memory region is usable at n samples.  
*/
static int good_range(char *start,memRange_t len,int n) {
  int i;
  memRange_t quant=divide_range(len,n);
  for (i=0;i<n;i++)
    if (bad_location(start+i*quant)) return 0; /* it's got some bad parts */
  /* It's all good: */
  return 1;
}

/*Check if this entire memory range, or some subset 
  of the range, is usable.  If so, write it into max.
*/
static void check_range(char *start,char *end,memRegion_t *max)
{
  memRange_t len;
  char *initialStart=start, *initialEnd=end;
  CmiUInt8 tb = (CmiUInt8)gig*1024ul;   /* One terabyte */
  CmiUInt8 vm_limit = tb*256ul;   /* terabyte */

  if (start>=end) return; /*Ran out of hole*/
  len=(memRange_t)end-(memRange_t)start;
  
#if 0
    /* too conservative */
  if (len/gig>64u) { /* This is an absurd amount of space-- cut it down, for safety */
     start+=16u*gig;
     end=start+32u*gig;
     len=(memRange_t)end-(memRange_t)start;  
  }
#else
  /* Note: 256TB == 248 bytes.  So a 48-bit virtual-address CPU 
   *    can only actually address 256TB of space. */
  if (len/tb>10u) { /* This is an absurd amount of space-- cut it down, for safety */
    const memRange_t other_libs=16ul*gig; /* space for other libraries to use */
    start+=other_libs;
    end=pmin(start+vm_limit-2*other_libs, end-other_libs);
    len=(memRange_t)end-(memRange_t)start;
  }
#endif
  if (len<=max->len) return; /*It's too short already!*/
#if CMK_THREADS_DEBUG
  CmiPrintf("[%d] Checking at %p - %p\n",CmiMyPe(),start,end);
#endif
  
  /* Check the middle of the range */
  if (!good_range(start,len,256)) {
    /* Try to split into subranges: */
    int i,n=2;
#if CMK_THREADS_DEBUG
    CmiPrintf("[%d] Trying to split bad address space at %p - %p...\n",CmiMyPe(),start,end);
#endif
    len=divide_range(len,n);
    for (i=0;i<n;i++) {
        char *cur=start+i*len;
	if (partially_good(cur,len,16))
	   check_range(cur,cur+len,max);
    }
    return; /* Hopefully one of the subranges will be any good */
  }
  else /* range is good */
  { 
#if CMK_THREADS_DEBUG
    CmiPrintf("[%d] Address space at %p - %p is largest\n",CmiMyPe(),start,end);
#endif

    /*If we got here, we're the new largest usable range*/
    max->len=len;
    max->start=start;
    max->type="Unused";
  }
}

/*Find the first available memory region of at least the
  given size not touching any data in the used list.
 */
static memRegion_t find_free_region(memRegion_t *used,int nUsed,int atLeast) 
{
  memRegion_t max;
  int i,j;  

  max.start=0; 
  max.len=atLeast;
  /*Find the largest hole between regions*/
  for (i=0;i<nUsed;i++) {
    /*Consider a hole starting at the end of region i*/
    char *holeStart=used[i].start+used[i].len;
    char *holeEnd=(void *)(-1);
    
    /*Shrink the hole by all others*/ 
    for (j=0;j<nUsed && pointer_lt(holeStart,holeEnd);j++) {
      if (pointer_lt(used[j].start,holeStart)) 
        holeStart=pmax(holeStart,used[j].start+used[j].len);
      else if (pointer_lt(used[j].start,holeEnd)) 
        holeEnd=pmin(holeEnd,used[j].start);
    } 

    check_range(holeStart,holeEnd,&max);
  }

  return max; 
}

/*
By looking at the address range carefully, try to find 
the largest usable free region on the machine.
*/
static int find_largest_free_region(memRegion_t *destRegion) {
    char *staticData =(char *) __static_data_loc();
    char *code = (char *)&find_free_region;
    char *threadData = (char *)&errno;
    char *codeDll = (char *)fprintf;
    char *heapLil = (char*) malloc(1);
    char *heapBig = (char*) malloc(6*meg);
    char *stack = (char *)__cur_stack_frame();
    size_t mmapAnyLen = 1*meg;
    char *mmapAny = (char*) call_mmap_anywhere(mmapAnyLen);

    int i,nRegions=8;
    memRegion_t regions[10]; /*used portions of address space*/
    memRegion_t freeRegion; /*Largest unused block of address space*/

/*Mark off regions of virtual address space as ususable*/
    regions[0].type="NULL";
    regions[0].start=NULL; regions[0].len=16u*meg;
    
    regions[1].type="Static program data";
    regions[1].start=staticData; regions[1].len=256u*meg;
    
    regions[2].type="Program executable code";
    regions[2].start=code; regions[2].len=256u*meg;
    
    regions[3].type="Heap (small blocks)";
    regions[3].start=heapLil; regions[3].len=1u*gig;
    
    regions[4].type="Heap (large blocks)";
    regions[4].start=heapBig; regions[4].len=1u*gig;
    
    regions[5].type="Stack space";
    regions[5].start=stack; regions[5].len=256u*meg;

    regions[6].type="Program dynamically linked code";
    regions[6].start=codeDll; regions[6].len=256u*meg; 

    regions[7].type="Result of a non-fixed call to mmap";
    regions[7].start=mmapAny; regions[7].len=1u*gig; 

    regions[8].type="Thread private data";
    regions[8].start=threadData; regions[8].len=256u*meg; 

    _MEMCHECK(heapBig); free(heapBig);
    _MEMCHECK(heapLil); free(heapLil); 
    call_munmap(mmapAny,mmapAnyLen);
    
    /*Align each memory region*/
    for (i=0;i<nRegions;i++) {
      memRegion_t old=regions[i];
      memRange_t p=(memRange_t)regions[i].start;
      p&=~(regions[i].len-1); /*Round start down to a len-boundary (mask off low bits)*/
      regions[i].start=(char *)p;
#if CMK_THREADS_DEBUG
      CmiPrintf("[%d] Memory map: %p - %p %s (at %p)\n",CmiMyPe(),
	      regions[i].start,regions[i].start+regions[i].len,
	      regions[i].type, old.start);
#endif
    }
    
    /*Find a large, unused region in this map: */
    freeRegion=find_free_region(regions,nRegions,(512u)*meg);
    
    if (freeRegion.start==0) 
    { /*No free address space-- disable isomalloc:*/
      return 0;
    }
    else /* freeRegion is valid */
    {
      *destRegion=freeRegion;
      
      return 1;
    }
}

static void init_ranges(char **argv)
{
  /*Largest value a signed int can hold*/
  memRange_t intMax=(((memRange_t)1)<<(sizeof(int)*8-1))-1;

  /*Round slot size up to nearest page size*/
  slotsize=16*1024;
  slotsize=(slotsize+CMK_MEMORY_PAGESIZE-1) & ~(CMK_MEMORY_PAGESIZE-1);
#if CMK_THREADS_DEBUG
  CmiPrintf("[%d] Using slotsize of %d\n", CmiMyPe(), slotsize);
#endif

  if (CmiMyRank()==0 && numslots==0)
  { /* Find the largest unused region of virtual address space */
      memRegion_t freeRegion;
      freeRegion.len=0u;
#ifdef CMK_MMAP_START_ADDRESS /* Hardcoded start address, for machines where automatic fails */
      freeRegion.start=CMK_MMAP_START_ADDRESS;
      freeRegion.len=CMK_MMAP_LENGTH_MEGS*meg;
#endif
      if (freeRegion.len==0u) find_largest_free_region(&freeRegion);
      
#if 0
      /*Make sure our largest slot number doesn't overflow an int:*/
      if (freeRegion.len/slotsize>intMax)
        freeRegion.len=intMax*slotsize;
#endif
      
      if (freeRegion.len==0u) {
        disable_isomalloc("no free virtual address space");
      }
      else /* freeRegion.len>0, so can isomalloc */
      {
#if CMK_THREADS_DEBUG
        CmiPrintf("[%d] Isomalloc memory region: %p - %p (%d megs)\n",CmiMyPe(),
	      freeRegion.start,freeRegion.start+freeRegion.len,
	      freeRegion.len/meg);
#endif
        
        /*Isomalloc covers entire unused region*/
        isomallocStart=freeRegion.start;
        isomallocEnd=freeRegion.start+freeRegion.len;
        numslots=(freeRegion.len/slotsize)/CmiNumPes();
        
#if CMK_THREADS_DEBUG
        CmiPrintf("[%d] Can isomalloc up to %lu megs per pe\n",CmiMyPe(),
	      ((memRange_t)numslots)*slotsize/meg);
#endif
      }
  }
  /*SMP Mode: wait here for rank 0 to initialize numslots so we can set up myss*/
  CmiNodeAllBarrier(); 
  
  if (isomallocStart!=NULL) {
    CpvInitialize(slotset *, myss);
    CpvAccess(myss) = new_slotset(pe2slot(CmiMyPe()), numslots);
  }
}


/************* Communication: for grabbing/freeing remote slots *********/
typedef struct _slotmsg
{
  char cmicore[CmiMsgHeaderSizeBytes];
  int pe; /*Source processor*/
  CmiInt8 slot; /*First requested slot*/
  CmiInt8 nslots; /*Number of requested slots*/
} slotmsg;

static slotmsg *prepare_slotmsg(CmiInt8 slot,CmiInt8 nslots)
{
	slotmsg *m=(slotmsg *)CmiAlloc(sizeof(slotmsg));
	m->pe=CmiMyPe();
	m->slot=slot;
	m->nslots=nslots;
	return m;
}

static void grab_remote(slotmsg *msg)
{
	grab_slots(CpvAccess(myss),msg->slot,msg->nslots);
	CmiFree(msg);
}

static void free_remote(slotmsg *msg)
{
	free_slots(CpvAccess(myss),msg->slot,msg->nslots);
	CmiFree(msg);
}
static int grab_remote_idx, free_remote_idx;

struct slotOP {
	/*Function pointer to perform local operation*/
	void (*local)(slotset *ss,CmiInt8 s,CmiInt8 n);
	/*Index to perform remote operation*/
	int remote;
};
typedef struct slotOP slotOP;
static slotOP grabOP,freeOP;

static void init_comm(char **argv)
{
	grab_remote_idx=CmiRegisterHandler((CmiHandler)grab_remote);
	free_remote_idx=CmiRegisterHandler((CmiHandler)free_remote);	
	grabOP.local=grab_slots;
	grabOP.remote=grab_remote_idx;
	freeOP.local=free_slots;
	freeOP.remote=free_remote_idx;	
}

/*Apply the given operation to the given slots which
  lie on the given processor.*/
static void one_slotOP(const slotOP *op,int pe,CmiInt8 s,CmiInt8 n)
{
/*Shrink range to only those covered by this processor*/
	/*First and last slot for this processor*/
	CmiInt8 p_s=pe2slot(pe), p_e=pe2slot(pe+1);
	CmiInt8 e=s+n;
	if (s<p_s) s=p_s;
	if (e>p_e) e=p_e;
	n=e-s;

/*Send off range*/
	if (pe==CmiMyPe()) 
		op->local(CpvAccess(myss),s,n);
	else 
	{/*Remote request*/
		slotmsg *m=prepare_slotmsg(s,n);
		CmiSetHandler(m, freeOP.remote);
		CmiSyncSendAndFree(pe,sizeof(slotmsg),m);
	}
}

/*Apply the given operation to all slots in the range [s, s+n) 
After a restart from checkpoint, a slotset can cross an 
arbitrary set of processors.
*/
static void all_slotOP(const slotOP *op,CmiInt8 s,CmiInt8 n)
{
	int spe=slot2pe(s), epe=slot2pe(s+n-1);
	int pe;
	for (pe=spe; pe<=epe; pe++)
		one_slotOP(op,pe,s,n);
}

/************** External interface ***************/
void *CmiIsomalloc(int size)
{
	CmiInt8 s,n;
	CmiIsomallocBlock *blk;
	if (isomallocStart==NULL) return disabled_map(size);
	n=length2slots(size);
	/*Always satisfy mallocs with local slots:*/
	s=get_slots(CpvAccess(myss),n);
	if (s==-1) {
		CmiError("Not enough address space left on processor %d to isomalloc %d bytes!\n",
			 CmiMyPe(),size);
		CmiAbort("Out of virtual address space for isomalloc");
	}
	grab_slots(CpvAccess(myss),s,n);
	blk=map_slots(s,n);
	if (!blk) map_failed(s,n);
	blk->slot=s;
	blk->length=size;
	return block2pointer(blk);
}

void CmiIsomallocPup(pup_er p,void **blockPtrPtr)
{
	CmiIsomallocBlock *blk;
	CmiInt8 s,length;
        CmiInt8 n;
	if (isomallocStart==NULL) CmiAbort("isomalloc is disabled-- cannot use IsomallocPup");

	if (!pup_isUnpacking(p)) 
	{ /*We have an existing block-- unpack start slot & length*/
		blk=pointer2block(*blockPtrPtr);
		s=blk->slot;
		length=blk->length;
	}
	
	pup_long(p,&s);
	pup_long(p,&length);
	n=length2slots(length);
	
	if (pup_isUnpacking(p)) 
	{ /*Must allocate a new block in its old location*/
		if (pup_isUserlevel(p) || pup_isRestarting(p))
			/*Checkpoint: must grab old slots (even remote!)*/
			all_slotOP(&grabOP,s,n);
		blk=map_slots(s,n);
		if (!blk) map_failed(s,n);
		blk->slot=s;
		blk->length=length;
		*blockPtrPtr=block2pointer(blk);
	}
	
	/*Pup the allocated data*/
	pup_bytes(p,*blockPtrPtr,length);
	
	if (pup_isDeleting(p)) 
	{ /*Unmap old slots, but do not mark as free*/
		unmap_slots(s,n);
		*blockPtrPtr=NULL; /*Zero out user's pointer*/
	}
}

void CmiIsomallocFree(void *blockPtr)
{
	if (isomallocStart==NULL) {
		disabled_unmap(blockPtr);
	}
	else if (blockPtr!=NULL)
	{
		CmiIsomallocBlock *blk=pointer2block(blockPtr);
		CmiInt8 s=blk->slot; 
                CmiInt8 n=length2slots(blk->length);
		unmap_slots(s,n);
		/*Mark used slots as free*/
		all_slotOP(&freeOP,s,n);
	}
}

CmiInt8   CmiIsomallocLength(void *block)
{
	return pointer2block(block)->length;
}

/*Return true if this address is in the region managed by isomalloc*/
int CmiIsomallocInRange(void *addr)
{
	if (isomallocStart==NULL) return 0; /* There is no range we manage! */
	return pointer_ge((char *)addr,isomallocStart) && 
	       pointer_lt((char*)addr,isomallocEnd);
}

void CmiIsomallocInit(char **argv)
{
#if CMK_NO_ISO_MALLOC
  disable_isomalloc("isomalloc disabled by conv-mach");
#else
  init_comm(argv);
  if (!init_map(argv)) {
    disable_isomalloc("mmap() does not work");
  }
  else {
    if (read_randomflag() == 1) {    /* randomization stack pointer */
      if (CmiMyPe() == 0)
        printf("Charm warning> Randomization of stack pointer is turned on in Kernel, run 'echo 0 > /proc/sys/kernel/randomize_va_space' as root to disable it. Thread migration may not work! \n");
    }
    init_ranges(argv);
  }
#endif
}

/***************** BlockList interface *********
This was moved here from memory-isomalloc.c when it 
was realized that a list-of-isomalloc'd-blocks is useful for
more than just isomalloc heaps.
*/

typedef CmiIsomallocBlockList Slot;

/*Convert a slot to a user address*/
static char *Slot_toUser(Slot *s) {return (char *)(s+1);}
static Slot *Slot_fmUser(void *s) {return ((Slot *)s)-1;}


/*Build a new blockList.*/
CmiIsomallocBlockList *CmiIsomallocBlockListNew(void)
{
	CmiIsomallocBlockList *ret;
	ret=(CmiIsomallocBlockList *)CmiIsomalloc(sizeof(*ret));
	ret->next=ret; /*1-entry circular linked list*/
	ret->prev=ret;
	return ret;
}

/*Pup all the blocks in this list.  This amounts to two circular
list traversals.  Because everything's isomalloc'd, we don't even
have to restore the pointers-- they'll be restored automatically!
*/
void CmiIsomallocBlockListPup(pup_er p,CmiIsomallocBlockList **lp)
{
	int i,nBlocks=0;
	Slot *cur=NULL, *start=*lp;
#if 0 /*#ifndef CMK_OPTIMIZE*/
	if (CpvAccess(isomalloc_blocklist)!=NULL)
		CmiAbort("Called CmiIsomallocBlockListPup while a blockList is active!\n"
			"You should swap out the active blocklist before pupping.\n");
#endif
	/*Count the number of blocks in the list*/
	if (!pup_isUnpacking(p)) {
		nBlocks=1; /*<- Since we have to skip the start block*/
		for (cur=start->next; cur!=start; cur=cur->next) 
			nBlocks++;
		/*Prepare for next trip around list:*/
		cur=start;
	}
	pup_int(p,&nBlocks);
	
	/*Pup each block in the list*/
	for (i=0;i<nBlocks;i++) {
		void *newBlock=cur;
		if (!pup_isUnpacking(p)) 
		{ /*While packing, we traverse the list to find our blocks*/
			cur=cur->next;
		}
		CmiIsomallocPup(p,&newBlock);
		if (i==0 && pup_isUnpacking(p))
			*lp=(Slot *)newBlock;
	}
	if (pup_isDeleting(p))
		*lp=NULL;
}

/*Delete all the blocks in this list.*/
void CmiIsomallocBlockListDelete(CmiIsomallocBlockList *l)
{
	Slot *start=l;
	Slot *cur=start;
	if (cur==NULL) return; /*Already deleted*/
	do {
		Slot *doomed=cur;
		cur=cur->next; /*Have to stash next before deleting cur*/
		CmiIsomallocFree(doomed);
	} while (cur!=start);
}

/*Allocate a block from this blockList*/
void *CmiIsomallocBlockListMalloc(CmiIsomallocBlockList *l,int nBytes)
{
	Slot *n; /*Newly created slot*/
	n=(Slot *)CmiIsomalloc(sizeof(Slot)+nBytes);
	/*Link the new block into the circular blocklist*/
	n->prev=l;
	n->next=l->next;
	l->next->prev=n;
	l->next=n;
	return Slot_toUser(n);
}

/*Remove this block from its list and memory*/
void CmiIsomallocBlockListFree(void *block)
{
	Slot *n=Slot_fmUser(block);
#if DOHEAPCHECK
	if (n->prev->next!=n || n->next->prev!=n) 
		CmiAbort("Heap corruption detected in isomalloc block list header!\n"
			"  Run with ++debug and look for writes to negative array indices");
#endif
	/*Link ourselves out of the blocklist*/
	n->prev->next=n->next;
	n->next->prev=n->prev;
	CmiIsomallocFree(n);
}





