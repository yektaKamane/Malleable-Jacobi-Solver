#define CONVERSE_VERSION_VMI                               1

#define CMK_CCS_AVAILABLE                                  1
#define NODE_0_IS_CONVHOST                                 1
#define CMK_WEB_MODE                                       0

#define CMK_CMIPRINTF_IS_A_BUILTIN                         0 

#define CMK_CMIDELIVERS_USE_COMMON_CODE                    1

#define CMK_PERSISTENT_COMM                                0

#if CMK_PERSISTENT_COMM
#define CMK_MULTICAST_GROUP_TYPE                struct { unsigned pe, id; }
#define CMK_MULTICAST_DEF_USE_COMMON_CODE                  1
#define CMK_MULTICAST_LIST_USE_COMMON_CODE                 0
#define CMK_MULTICAST_LIST_USE_SPECIAL_CODE                1
#define CMK_MULTICAST_GROUP_USE_COMMON_CODE                1
#else
#define CMK_MULTICAST_GROUP_TYPE                struct { unsigned pe, id; }
#define CMK_MULTICAST_DEF_USE_COMMON_CODE                  1
#define CMK_MULTICAST_LIST_USE_COMMON_CODE                 1
#define CMK_MULTICAST_LIST_USE_SPECIAL_CODE                0
#define CMK_MULTICAST_GROUP_USE_COMMON_CODE                1
#endif

#define CMK_BROADCAST_SPANNING_TREE                        1
#define CMK_SPANTREE_MAXSPAN                               4
#define CMK_SPANTREE_USE_COMMON_CODE                       1

#define CMK_MSG_HEADER_BASIC  CMK_MSG_HEADER_EXT
#define CMK_MSG_HEADER_EXT    { CmiUInt2 vmitype,tree_rank,tree_root,hdl,xhdl,info,stratid,root; }
#define CMK_MSG_HEADER_BLUEGENE    { CmiUInt2 d0,d1,d2,d3,d4,d5,hdl,xhdl,info,stratid,root,padding3; int nd, n; double rt; CmiInt2 tID; CmiUInt2 hID; char t; int msgID; int srcPe;}

#define CMK_VECTOR_SEND_USES_COMMON_CODE                   1

#define CMK_IMMEDIATE_MSG                                  0

#define CMK_ASYNC_NOT_NEEDED                               0
#define CMK_ASYNC_USE_FIOASYNC_AND_FIOSETOWN               0
#define CMK_ASYNC_USE_FIOASYNC_AND_SIOCSPGRP               0
#define CMK_ASYNC_USE_FIOSSAIOSTAT_AND_FIOSSAIOOWN         0
#define CMK_ASYNC_USE_F_SETFL_AND_F_SETOWN                 1

#define CMK_GETPAGESIZE_AVAILABLE                          1

#define CMK_MALLOC_USE_GNU_MALLOC                          0
#define CMK_MALLOC_USE_OS_BUILTIN                          1
#define CMK_MALLOC_USE_GNUOLD_MALLOC                       0

#define CMK_MEMORY_PROTECTABLE                             1

#define CMK_NODE_QUEUE_AVAILABLE                           0

#define CMK_RSH_IS_A_COMMAND                               1
#define CMK_RSH_NOT_NEEDED                                 0
#define CMK_RSH_USE_REMSH                                  0

#define CMK_SHARED_VARS_EXEMPLAR                           0
#define CMK_SHARED_VARS_UNAVAILABLE                        1
#define CMK_SHARED_VARS_SUN_THREADS                        0
#define CMK_SHARED_VARS_UNIPROCESSOR                       0

#define CMK_THREADS_COPY_STACK                             0
#define CMK_THREADS_USE_CONTEXT                            1
#define CMK_THREADS_USE_PTHREADS                           0
#define CMK_THREADS_ARE_WIN32_FIBERS                       0

#define CMK_THREADS_REQUIRE_NO_CPV                         0

#define CMK_SIGNAL_NOT_NEEDED                              0
#define CMK_SIGNAL_USE_SIGACTION                           1
#define CMK_SIGNAL_USE_SIGACTION_WITH_RESTART              0

#define CMK_TIMER_USE_RDTSC                                1
#define CMK_TIMER_USE_GETRUSAGE                            0
#define CMK_TIMER_USE_SPECIAL                              0
#define CMK_TIMER_USE_TIMES                                0

#define CMK_TYPEDEF_INT2 short
#define CMK_TYPEDEF_INT4 int
#define CMK_TYPEDEF_INT8 long long
#define CMK_TYPEDEF_UINT2 unsigned short
#define CMK_TYPEDEF_UINT4 unsigned int
#define CMK_TYPEDEF_UINT8 unsigned long long
#define CMK_TYPEDEF_FLOAT4 float
#define CMK_TYPEDEF_FLOAT8 double

#define CMK_WHEN_PROCESSOR_IDLE_BUSYWAIT                   1
#define CMK_WHEN_PROCESSOR_IDLE_USLEEP                     0

#define CMK_LBDB_ON                                        1
