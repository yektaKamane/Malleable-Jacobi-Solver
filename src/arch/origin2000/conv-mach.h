/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef _CONV_MACH_H
#define _CONV_MACH_H


#define CMK_ORIGIN2000					   1
#define CMK_CCS_AVAILABLE                                  1

#define CMK_CMIDELIVERS_USE_COMMON_CODE                    1
#define CMK_CMIDELIVERS_USE_SPECIAL_CODE                   0

#define CMK_CMIPRINTF_IS_A_BUILTIN                         0
#define CMK_CMIPRINTF_IS_JUST_PRINTF                       1

#define CMK_DEFAULT_MAIN_USES_COMMON_CODE                  1


#define CMK_GETPAGESIZE_AVAILABLE                          1

#define CMK_HANDLE_SIGUSR                                  1

#define CMK_IS_HETERO                                      0


#define CMK_MALLOC_USE_GNU_MALLOC                          0
#define CMK_MALLOC_USE_OS_BUILTIN                          1

#define CMK_MEMORY_PAGESIZE                                8192
#define CMK_MEMORY_PROTECTABLE                             0

#define CMK_MSG_HEADER_BASIC  CMK_MSG_HEADER_EXT
#define CMK_MSG_HEADER_EXT    { CmiUInt2 hdl,xhdl,info,padding; }
#define CMK_MSG_HEADER_BLUEGENE    { CmiUInt2 hdl,xhdl,info,padding; int nd, n; double rt; CmiInt2 tID; CmiUInt2 hID; char t; int msgID; int srcPe;}

#define CMK_MULTICAST_GROUP_TYPE                struct { unsigned pe, id; }
#define CMK_MULTICAST_DEF_USE_COMMON_CODE                  1
#define CMK_MULTICAST_LIST_USE_COMMON_CODE                 1
#define CMK_MULTICAST_GROUP_USE_COMMON_CODE                1

#define CMK_NODE_QUEUE_AVAILABLE                           0

#define CMK_REDUCTION_USES_COMMON_CODE                     1
#define CMK_REDUCTION_USES_SPECIAL_CODE                    0

#define CMK_RSH_IS_A_COMMAND                               0
#define CMK_RSH_NOT_NEEDED                                 1
#define CMK_RSH_USE_REMSH                                  0

#define CMK_SHARED_VARS_EXEMPLAR                           0
#define CMK_SHARED_VARS_UNAVAILABLE                        1
#define CMK_SHARED_VARS_UNIPROCESSOR                       0
#define CMK_SHARED_VARS_SUN_THREADS                        0

#define CMK_THREADS_USE_CONTEXT                            0
#define CMK_CONTEXT_STACKEND                               1

#define CMK_SIGNAL_NOT_NEEDED                              0
#define CMK_SIGNAL_USE_SIGACTION                           0
#define CMK_SIGNAL_USE_SIGACTION_WITH_RESTART              1

#define CMK_SPANTREE_MAXSPAN                               4
#define CMK_SPANTREE_USE_COMMON_CODE                       1
#define CMK_SPANTREE_USE_SPECIAL_CODE                      0



#define CMK_SYNCHRONIZE_ON_TCP_CLOSE                       0

#define CMK_THREADS_REQUIRE_NO_CPV                         0
#define CMK_THREADS_COPY_STACK                             0

#define CMK_TIMER_USE_GETRUSAGE                            0
#define CMK_TIMER_USE_SPECIAL                              0
#define CMK_TIMER_USE_TIMES                                0

#define CMK_TYPEDEF_INT2 short
#define CMK_TYPEDEF_INT4 int
#define CMK_TYPEDEF_INT8 long
#define CMK_TYPEDEF_UINT2 unsigned short
#define CMK_TYPEDEF_UINT4 unsigned int
#define CMK_TYPEDEF_UINT8 unsigned long
#define CMK_TYPEDEF_FLOAT4 float
#define CMK_TYPEDEF_FLOAT8 double


#define CMK_VECTOR_SEND_USES_COMMON_CODE                   1
#define CMK_VECTOR_SEND_USES_SPECIAL_CODE                  0


#define CMK_WHEN_PROCESSOR_IDLE_BUSYWAIT                   1
#define CMK_WHEN_PROCESSOR_IDLE_USLEEP                     0


#define NODE_0_IS_CONVHOST                                 1
#define CMK_DEBUG_MODE                                     0
#define CMK_WEB_MODE                                       1


#define CMK_LBDB_ON					   1




#endif

