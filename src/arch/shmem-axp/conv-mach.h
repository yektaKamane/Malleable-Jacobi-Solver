/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef _CONV_MACH_H
#define _CONV_MACH_H

#define CMK_CCS_AVAILABLE                                  1

#define CMK_CMIDELIVERS_USE_COMMON_CODE                    1

#define CMK_CMIPRINTF_IS_A_BUILTIN                         0

#define CMK_GETPAGESIZE_AVAILABLE                          0

#define CMK_HANDLE_SIGUSR                                  1

#define CMK_IS_HETERO                                      0


#define CMK_MALLOC_USE_GNU_MALLOC                          0
#define CMK_MALLOC_USE_OS_BUILTIN                          1

#define CMK_MEMORY_PAGESIZE                                8192
#define CMK_MEMORY_PROTECTABLE                             0

#define CMK_MSG_HEADER_BASIC  CMK_MSG_HEADER_EXT
#define CMK_MSG_HEADER_EXT_   char gap[56]; CmiUInt2 hdl,xhdl,info,stratid,root,padding1,padding2,padding3;
#define CMK_MSG_HEADER_EXT       { CMK_MSG_HEADER_EXT_ }
#define CMK_MSG_HEADER_BLUEGENE  {CMK_MSG_HEADER_EXT_ CMK_BLUEGENE_FIELDS}

#define CMK_MULTICAST_GROUP_TYPE                struct { unsigned pe, id; }
#define CMK_MULTICAST_DEF_USE_COMMON_CODE                  1
#define CMK_MULTICAST_LIST_USE_COMMON_CODE                 0
#define CMK_MULTICAST_GROUP_USE_COMMON_CODE                0

#define CMK_NODE_QUEUE_AVAILABLE                           0

#define CMK_RSH_IS_A_COMMAND                               0
#define CMK_RSH_NOT_NEEDED                                 1
#define CMK_RSH_USE_REMSH                                  0

#define CMK_SHARED_VARS_EXEMPLAR                           0
#define CMK_SHARED_VARS_UNAVAILABLE                        1
#define CMK_SHARED_VARS_UNIPROCESSOR                       0

#define CMK_SIGNAL_NOT_NEEDED                              0
#define CMK_SIGNAL_USE_SIGACTION                           0
#define CMK_SIGNAL_USE_SIGACTION_WITH_RESTART              1

#define CMK_SPANTREE_MAXSPAN                               4
#define CMK_SPANTREE_USE_COMMON_CODE                       1

#define CMK_SYNCHRONIZE_ON_TCP_CLOSE                       0

#define CMK_THREADS_REQUIRE_NO_CPV                         0
#define CMK_THREADS_COPY_STACK                             0

#define CMK_TIMER_USE_GETRUSAGE                            1
#define CMK_TIMER_USE_SPECIAL                              0
#define CMK_TIMER_USE_TIMES                                0

#define CMK_TYPEDEF_INT2 short
#define CMK_TYPEDEF_INT4 short
#define CMK_TYPEDEF_INT8 int
#define CMK_TYPEDEF_UINT2 unsigned short
#define CMK_TYPEDEF_UINT4 unsigned short
#define CMK_TYPEDEF_UINT8 unsigned int
#define CMK_TYPEDEF_FLOAT4 float
#define CMK_TYPEDEF_FLOAT8 double

#define CMK_VECTOR_SEND_USES_COMMON_CODE                        1

#define CMK_WHEN_PROCESSOR_IDLE_BUSYWAIT                   1
#define CMK_WHEN_PROCESSOR_IDLE_USLEEP                     0

#define CMK_USE_HP_MAIN_FIX                                1

#define NODE_0_IS_CONVHOST                                 1
#define CMK_DEBUG_MODE                                     0
#define CMK_WEB_MODE                                       1


#define CMK_LBDB_ON					   1

#define CMK_SHMEM_H					<shmem.h>
#define CMK_SHMEM_INIT					   shmem_init()
#define CMK_SHMEM_LOCK					   1

#define CMK_TRACE_LOGFILE_NUM_CONTROL                      0


#endif

