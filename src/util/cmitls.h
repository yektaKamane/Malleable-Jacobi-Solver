#if !defined(CMITLS_H)
#define CMITLS_H

#include "conv-config.h"

#if CMK_HAS_TLS_VARIABLES

#if CMK_HAS_ELF_H \
    && ((CMK_DLL_USE_DLOPEN && CMK_HAS_RTLD_DEFAULT) || CMK_HAS_DL_ITERATE_PHDR)

#include <elf.h>

#if ( defined(__LP64__) || defined(_LP64) )
#define ELF64
#else
#define ELF32
#endif

#endif

#endif

#if defined ELF32
typedef Elf32_Addr Addr;
typedef Elf32_Ehdr Ehdr;
typedef Elf32_Phdr Phdr;
#elif defined ELF64
typedef Elf64_Addr Addr;
typedef Elf64_Ehdr Ehdr;
typedef Elf64_Phdr Phdr;
#else
typedef void * Addr;
#endif

typedef struct tlsseg {
  Addr memseg;
  size_t size;
  size_t align;
} tlsseg_t;

#ifdef __cplusplus
extern "C" {
#endif

void CmiTLSInit(void);
void allocNewTLSSeg(tlsseg_t* t, CthThread th);
void switchTLS(tlsseg_t*, tlsseg_t*);
void currentTLS(tlsseg_t*);

#ifdef __cplusplus
}
#endif

#endif
