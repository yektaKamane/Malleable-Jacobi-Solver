/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef QUEUEING_H
#define QUEUEING_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CINTBITS
#define CINTBITS ((unsigned int) (sizeof(int)*8))
#endif
typedef struct prio_struct
{
  unsigned short bits;
  unsigned short ints;
  unsigned int data[1];
}
*prio;

typedef struct deq_struct
{
  /* Note: if head==tail, circ is empty */
  void **bgn; /* Pointer to first slot in circular buffer */
  void **end; /* Pointer past last slot in circular buffer */
  void **head; /* Pointer to first used slot in circular buffer */
  void **tail; /* Pointer to next available slot in circular buffer */
  void *space[4]; /* Enough space for the first 4 entries */
}
*deq;

typedef struct prioqelt_struct
{
  struct deq_struct data;
  struct prioqelt_struct *ht_next; /* Pointer to next bucket in hash table. */
  struct prioqelt_struct **ht_handle; /* Pointer to pointer that points to me (!) */
  struct prio_struct pri;
}
*prioqelt;

#define PRIOQ_TABSIZE 1017
typedef struct prioq_struct
{
  int heapsize;
  int heapnext;
  prioqelt *heap;
  prioqelt hashtab[PRIOQ_TABSIZE];
}
*prioq;

typedef struct Queue_struct
{
  unsigned int length;
  unsigned int maxlen;
  struct deq_struct zeroprio;
  struct prioq_struct negprioq;
  struct prioq_struct posprioq;
}
*Queue;

Queue CqsCreate(void);
void CqsDelete(Queue);
void CqsEnqueue(Queue, void *msg);
void CqsEnqueueFifo(Queue, void *msg);
void CqsEnqueueLifo(Queue, void *msg);
void CqsEnqueueGeneral(Queue, void *msg,int strategy, 
	       int priobits, unsigned int *prioPtr);

void CqsEnumerateQueue(Queue q, void ***resp);
void CqsDequeue(Queue, void **msgPtr);

unsigned int CqsLength(Queue);
unsigned int CqsMaxLength(Queue);
int CqsEmpty(Queue);
int  CqsPrioGT(prio, prio);
prio CqsGetPriority(Queue);

#ifdef __cplusplus
};
#endif

#endif
