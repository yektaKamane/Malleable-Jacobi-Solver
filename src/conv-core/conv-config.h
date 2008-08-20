/*
  Generic converse configuration-flags header.
*/
#ifndef __UIUC_CHARM_CONV_CONFIG_H
#define __UIUC_CHARM_CONV_CONFIG_H

/* 
 Include the automatically determined options.
  conv-autoconfig.h is written by the configure script.

 This header declares all automatically-determined properties
 of the machine, like compiler features, headers, and syscalls.
 Generally, more features should be moved into this header
 from the other manually generated headers.
*/
#include "conv-autoconfig.h"

/* 
 Include the machine.c configuration header 
  (e.g., charm/src/arch/net/conv-common.h )

 This header declares communication properties
 (like message header formats) and the various 
 machine arcana handled by machine.c.
*/ 
#include "conv-common.h"

/* 
 Include the system/platform header.
  (e.g., charm/src/arch/net-linux/conv-mach.h )
 
 This header declares the handling of malloc (OS or charm),
 signals, threads, timers, and other details. 
*/
#include "conv-mach.h"

/* 
 Include the build-time options.
  conv-mach-opt.h is written by the build script.

 This header includes any special build-time options.
 It's typically empty or very short.
*/
#include "conv-mach-opt.h"


/* Fix various invariants (these symbols should probably 
  be eliminated entirely) */
#define CMK_LBDB_OFF (!CMK_LBDB_ON)

#ifndef CMK_USE_HP_MAIN_FIX
# define CMK_USE_HP_MAIN_FIX	0
#endif

#if CMK_AMD64 && !defined(CMK_64BIT)
#define CMK_64BIT		1
#endif

#if CMK_BLUEGENEL
#define CMK_VERSION_BLUEGENE	1
#endif

#endif
