
#define CMK_CMIDELIVERS_USE_COMMON_CODE                    1
#define CMK_CMIDELIVERS_USE_SPECIAL_CODE                   0

#define CMK_CMIPRINTF_IS_A_BUILTIN                         0
#define CMK_CMIPRINTF_IS_JUST_PRINTF                       1

#define CMK_HANDLE_SIGUSR                                  1

#define CMK_MSG_HEADER_BASIC  CMK_MSG_HEADER_EXT
#define CMK_MSG_HEADER_EXT    { CmiUInt2 rank,root,hdl,xhdl,info,stratid,padding1,padding2;}
#define CMK_MSG_HEADER_BLUEGENE    { CmiUInt2 rank,root,hdl,xhdl,info,stratid,padding1,padding2; int nd, n; double rt; CmiInt2 tID; CmiUInt2 hID; char t; int msgID; int srcPe;}

#define CMK_MULTICAST_GROUP_TYPE                struct { unsigned pe, id; }
#define CMK_MULTICAST_DEF_USE_COMMON_CODE                  1
#define CMK_MULTICAST_LIST_USE_COMMON_CODE                 1
#define CMK_MULTICAST_GROUP_USE_COMMON_CODE                1

#define CMK_REDUCTION_USES_COMMON_CODE                     1
#define CMK_REDUCTION_USES_SPECIAL_CODE                    0

#define CMK_RSH_IS_A_COMMAND                               0
#define CMK_RSH_NOT_NEEDED                                 1
#define CMK_RSH_USE_REMSH                                  0

#define CMK_SPANTREE_MAXSPAN                               4
#define CMK_SPANTREE_USE_COMMON_CODE                       1
#define CMK_SPANTREE_USE_SPECIAL_CODE                      0

#define CMK_VECTOR_SEND_USES_COMMON_CODE                   1
#define CMK_VECTOR_SEND_USES_SPECIAL_CODE                  0

#define CMK_CCS_AVAILABLE                                  1

#define NODE_0_IS_CONVHOST                                 1

#define CMK_IMMEDIATE_MSG				   0
