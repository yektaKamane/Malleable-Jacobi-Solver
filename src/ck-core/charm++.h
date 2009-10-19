/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef _CHARMPP_H_
#define _CHARMPP_H_

#include <stdlib.h>
#include <memory.h>
#include "charm.h"

#include "middle.h"

class CMessage_CkArgMsg {
public: static int __idx;
};
#define CK_ALIGN(val,to) (((val)+(to)-1)&~((to)-1))

#include "pup.h"
#include "cklists.h"
#include "ckbitvector.h"
#include "init.h"
#include "debug-charm.h"
#include "simd.h"

PUPbytes(CkChareID)
PUPbytes(CkGroupID)

  
/**
 * CkMessage is the superclass of all Charm++ messages.
 * Typically, a message foo inherits from CMessage_foo, which
 * inherits from CkMessage.  In the internals of Charm++,
 * messages are often represented by bare "void *"s, which is 
 * silly and dangerous.
 */
class CkMessage { 
	//Don't use these: use CkCopyMsg
	CkMessage(const CkMessage &);
	void operator=(const CkMessage &);
public:
	CkMessage() {}
	void operator delete(void *ptr) { CkFreeMsg(ptr); }
	
	/* This pup routine only packs the message itself, *not* the
	message header.  Use CkPupMessage instead of calling this directly. */
	void pup(PUP::er &p);
	
	/// This is used to display message contents in the debugger.
	static void ckDebugPup(PUP::er &p,void *msg);
};
class CMessage_CkMessage {
public:
	static int __idx;
};

/// CkArgMsg is passed to the mainchare's constructor.
class CkArgMsg : public CkMessage {
public:
  int argc;
  char **argv;
};

class CkArray;

void CkPupMessage(PUP::er &p,void **atMsg,int pack_detail=1);

//This is for passing a single Charm++ message via parameter marshalling
class CkMarshalledMessage {
	void *msg;
	//Don't use these: only pass by reference
	void operator=(const CkMarshalledMessage &);
 public:
	CkMarshalledMessage(void): msg(NULL) {}
	CkMarshalledMessage(CkMessage *m): msg(m) {} //Takes ownership of message
	CkMarshalledMessage(const CkMarshalledMessage &);
	~CkMarshalledMessage() {if (msg) CkFreeMsg(msg);}
	CkMessage *getMessage(void) {void *ret=msg; msg=NULL; return (CkMessage *)ret;}
	void pup(PUP::er &p) {CkPupMessage(p,&msg,1);}
};
PUPmarshall(CkMarshalledMessage)

/**
 * CkEntryOptions describes the options associated
 * with an entry method invocation, which include
 * the message priority and queuing strategy.
 * It is only used with parameter marshalling.
 */
class CkEntryOptions : public CkNoncopyable {
	int queueingtype; //CK_QUEUEING type
	int prioBits; //Number of bits of priority to use
	typedef unsigned int prio_t; //Datatype used to represent priorities
	prio_t *prioPtr; //Points to message priority values
	prio_t prioStore; //For short priorities, stores the priority value
public:
	CkEntryOptions(void): queueingtype(CK_QUEUEING_FIFO), prioBits(0), 
                              prioPtr(NULL), prioStore(0) {}

	~CkEntryOptions() {
		if ( prioPtr != NULL && queueingtype != CK_QUEUEING_IFIFO ) {
			delete [] prioPtr;
			prioBits = 0;
		}
	}
	
	inline void setPriority(prio_t integerPrio) {
		queueingtype=CK_QUEUEING_IFIFO;
		prioBits=8*sizeof(integerPrio);
		prioPtr=&prioStore;
		prioStore=integerPrio;
	}
	inline void setPriority(int prioBits_,const prio_t *prioPtr_) {
		if ( prioPtr != NULL && queueingtype != CK_QUEUEING_IFIFO ) {
			delete [] prioPtr;
			prioBits = 0;
		}
		queueingtype=CK_QUEUEING_BFIFO;
		prioBits=prioBits_;
		int dataLength = (prioBits + (sizeof(prio_t)*8 - 1)) /
		                 (sizeof(prio_t)*8);
		prioPtr = new prio_t[dataLength];
		memcpy((void *)prioPtr, prioPtr_, dataLength*sizeof(unsigned int));
	}
	inline void setPriority(const CkBitVector &cbv) {
		if ( cbv.data != NULL ) {
			if ( prioPtr != NULL && queueingtype != CK_QUEUEING_IFIFO ) {
				delete [] prioPtr;
				prioBits = 0;
			}
			queueingtype=CK_QUEUEING_BFIFO;
			prioBits=cbv.usedBits;
			int dataLength = (prioBits + (sizeof(prio_t)*8 - 1)) /
		                 	(sizeof(prio_t)*8);
			prioPtr = new prio_t[dataLength];
			memcpy((void *)prioPtr, cbv.data, dataLength*sizeof(prio_t));
		} else {
			queueingtype=CK_QUEUEING_BFIFO;
			prioBits=0;
			int dataLength = 1;
			prioPtr = new prio_t[dataLength];
			prioPtr[0] = 0;
		}
	}
	
	inline void setQueueing(int queueingtype_) {queueingtype=queueingtype_;}

	///These are used by CkAllocateMarshallMsg, below:
	inline int getQueueing(void) const {return queueingtype;}
	inline int getPriorityBits(void) const {return prioBits;}
	inline const prio_t *getPriorityPtr(void) const {return prioPtr;}
};

#include "CkMarshall.decl.h"
//This is the message type marshalled parameters get packed into:
class CkMarshallMsg : public CMessage_CkMarshallMsg {
public: 
	char *msgBuf;
};



//A queue-of-messages, like CkMsgQ<CkReductionMsg>
template <class MSG>
class CkMsgQ : public CkQ<MSG *> {
public:
	~CkMsgQ() { //Delete the messages in the queue:
		MSG *m;
		while (NULL!=(m=this->deq())) delete m;
	}
	void pup(PUP::er &p) {
		int l=this->length();
		p(l);
		for (int i=0;i<l;i++) {
			MSG *m=NULL;
			if (!p.isUnpacking()) m=this->deq();
			CkPupMessage(p,(void **)&m);
			this->enq(m);
		}
	}
	friend void operator|(PUP::er &p,CkMsgQ<MSG> &v) {v.pup(p);}
};

/*******************************************************
Array Index class.  An array index is just a hash key-- 
a run of integers used to look up an object in a hash table.
*/

#include "ckhashtable.h"

#ifndef CK_ARRAYINDEX_MAXLEN 
#define CK_ARRAYINDEX_MAXLEN 3 /*Max. # of integers in an array index*/
#endif

class CkArrayIndex
{
public:
        ///Length of index in *integers*
	short int nInts;
	///Number of dimensions in this index, not valid for user-defined indices
	short int dimension; 
	
	//Index data immediately follows...
	
	int *data(void) {return (int*)((&nInts)+2);}
	const int *data(void) const {return (int*)((&nInts)+2);}

        int getCombinedCount(void) const {
          if (nInts == 1) return data()[0];
          else if (nInts == 2) return data()[0] * data()[1];
          else if (nInts == 3) return data()[0] * data()[1] * data()[2];
          else return 0;
        }

	void pup(PUP::er &p);

    //These routines allow CkArrayIndex to be used in
    //  a CkHashtableT
	CkHashCode hash(void) const;
	static CkHashCode staticHash(const void *a,size_t);
	int compare(const CkArrayIndex &ind) const;
	static int staticCompare(const void *a,const void *b,size_t);
};

inline CkHashCode CkArrayIndex::hash(void) const
{
        register int i;
	register const int *d=data();
	register CkHashCode ret=d[0];
	for (i=1;i<nInts;i++)
		ret +=circleShift(d[i],10+11*i)+circleShift(d[i],9+7*i);
	return ret;
}
inline int CkArrayIndex::compare(const CkArrayIndex &i2) const
{
	const CkArrayIndex &i1=*this;
#if CMK_1D_ONLY
	return i1.data()[0]==i2.data()[0];
#else
	const int *d1=i1.data();
	const int *d2=i2.data();
	int l=i1.nInts;
	if (l!=i2.nInts) return 0;
	for (int i=0;i<l;i++)
		if (d1[i]!=d2[i])
			return 0;
	//If we got here, the two keys must have exactly the same data
	return 1;
#endif
}


//This class is as large as any CkArrayIndex
class CkArrayIndexMax : public CkArrayIndex {
	struct {
		int data[CK_ARRAYINDEX_MAXLEN];
	} index;
	void copyFrom(const CkArrayIndex &that)
	{
		nInts=that.nInts;
		dimension=that.dimension;
		index=((const CkArrayIndexMax *)&that)->index;
		//for (int i=0;i<nInts;i++) index.data[i]=that.data()[i];
	}
public:
	CkArrayIndexMax(void) { nInts=0; dimension=0; }
	CkArrayIndexMax(int i) { nInts=0; dimension=0; }   // used for CkVec
	CkArrayIndexMax(const CkArrayIndex &that) 
		{copyFrom(that);}
	CkArrayIndexMax &operator=(const CkArrayIndex &that) 
		{copyFrom(that); return *this;}
        void print() { CmiPrintf("%d: %d %d %d\n", nInts,index.data[0], index.data[1], index.data[2]); }
	void pup(PUP::er &p) {
		p|nInts;
		p|dimension;
		for (int i=0;i<nInts;i++) p|index.data[i];
	}
	/*
	 * Code for the previous attempt to introduce CkArrayID into CmiObjId
	 * 
	 * Turns out this:
	 * i) does not appear to work right.
	 * ii) does not result in the same ID as the ones used by the Load
	 *     Balancer.
	 *
	 * FIXME **CWL** - temporary solution:
	 * i) ignore the arrayID passed into the method.
	 * ii) use exactly the same code as the 
	 *     idx2LDObjid(const CkArrayIndex &idx) found in cklocation.C
	 *     The only problem is that the code only exists when the
	 *     CMK_LBDB_ON is set, so we cannot reuse the code and
	 *     risk having to modify the code should idx2LDObjid be changed.
	 *
	CmiObjId *getProjectionID(int arrayID) { 
	  CmiObjId *ret = new CmiObjId;
	  for (int i=0; i<CK_ARRAYINDEX_MAXLEN; i++) {
	    ret->id[i] = index.data[i];
	  }
	  ret->id[CK_ARRAYINDEX_MAXLEN] = arrayID;
	  return ret; 
	}
	*/
	CmiObjId *getProjectionID(int arrayID) {
	  CmiObjId *ret = new CmiObjId;
	    int i;
	    const int *data=this->data();
	    if (OBJ_ID_SZ>=this->nInts) {
	      for (i=0;i<this->nInts;i++)
		ret->id[i]=data[i];
	      for (i=this->nInts;i<OBJ_ID_SZ;i++)
		ret->id[i]=0;
	    } else {
	      //Must hash array index into LBObjid
	      int j;
	      for (j=0;j<OBJ_ID_SZ;j++)
		ret->id[j]=data[j];
	      for (i=0;i<this->nInts;i++)
		for (j=0;j<OBJ_ID_SZ;j++)
		  ret->id[j]+=circleShift(data[i],22+11*i*(j+1))+
		    circleShift(data[i],21-9*i*(j+1));
	    }
	    return ret;
	}

        CmiBool operator==(const CkArrayIndexMax& idx) const {
          if (nInts != idx.nInts) return CmiFalse;
          for (int i=0; i<nInts; i++)
                if (index.data[i] != idx.index.data[i]) return CmiFalse;
          return CmiTrue;
        }
};
PUPmarshall(CkArrayIndexMax)

//A layout-compatible version of a CkArrayIndexMax.
//  Needed, e.g., for use in unions where a constructor is forbidden.
class CkArrayIndexStruct {
public:
	short int nInts;
	short int dimension;
	int index[CK_ARRAYINDEX_MAXLEN];
	CkArrayIndexMax &asMax(void) 
		{return *(CkArrayIndexMax *)this;}
	const CkArrayIndexMax &asMax(void) const
		{return *(const CkArrayIndexMax *)this;}
	void pup(PUP::er &p) {
		p|nInts;
		p|dimension;
		for (int i=0;i<nInts;i++) p|index[i];
	}
};
PUPmarshall(CkArrayIndexStruct)

class CkArrayID {
	CkGroupID _gid;
public:
	CkArrayID() : _gid() { }
	CkArrayID(CkGroupID g) :_gid(g) {}
	inline void setZero(void) {_gid.setZero();}
	inline int isZero(void) const {return _gid.isZero();}
	operator CkGroupID() const {return _gid;}
	CkArray *ckLocalBranch(void) const
		{ return (CkArray *)CkLocalBranch(_gid); }
	static CkArray *CkLocalBranch(CkArrayID id) 
		{ return (CkArray *)::CkLocalBranch(id); }
	void pup(PUP::er &p) {p | _gid; }
	int operator == (const CkArrayID& other) const {
		return (_gid == other._gid);
	}
};
PUPmarshall(CkArrayID)

#include "cksection.h"

#include "ckcallback.h"

/********************* Superclass of all Chares ******************/
#if CMK_MULTIPLE_DELETE
#define CHARM_INPLACE_NEW \
    void *operator new(size_t, void *ptr) { return ptr; }; \
    void operator delete(void*, void*) {}; \
    void *operator new(size_t s) { return malloc(s); } \
    void operator delete(void *ptr) { free(ptr); }
#else
#define CHARM_INPLACE_NEW \
    void *operator new(size_t, void *ptr) { return ptr; }; \
    void *operator new(size_t s) { return malloc(s); } \
    void operator delete(void *ptr) { free(ptr); }
#endif

// for object message queue
#include "ckobjQ.h"


#ifdef _FAULT_MLOG_
class ChareMlogData;
#endif


/**
  The base class of all parallel objects in Charm++,
  including Array Elements, Groups, and NodeGroups.
*/
class Chare {
  protected:
    CkChareID thishandle;
#if CMK_OBJECT_QUEUE_AVAILABLE
    CkObjectMsgQ objQ;                // object message queue
#endif
  public:
#if CMK_FT_CHARE
    int chareIdx;                  // index in the chare obj table (chare_objs)
#endif
#ifdef _FAULT_MLOG_
    ChareMlogData *mlogData;
#endif
    Chare(CkMigrateMessage *m);
    Chare();
    virtual ~Chare(); //<- needed for *any* child to have a virtual destructor
    virtual void pup(PUP::er &p);//<- pack/unpack routine
    inline const CkChareID &ckGetChareID(void) const {return thishandle;}
    inline void CkGetChareID(CkChareID *dest) const {*dest=thishandle;}
    // object message queue
    void  CkEnableObjQ();
#if CMK_OBJECT_QUEUE_AVAILABLE
    inline CkObjectMsgQ &CkGetObjQueue() { return objQ; }
#endif
    CHARM_INPLACE_NEW
    /// Return the type of this chare, as present in _chareTable
    virtual int ckGetChareType() const;
    /// Return a strdup'd array containing this object's string name.
    virtual char *ckDebugChareName(void);
    /// Place into str a copy of the id of this object up to limit bytes, return
    /// the number of bytes used for the id
    virtual int ckDebugChareID(char *str, int limit);
    virtual void ckDebugPup(PUP::er &p);
    /// Called when a [threaded] charm entry method is created:
    virtual void CkAddThreadListeners(CthThread tid, void *msg);
};

//Superclass of all Groups that cannot participate in reductions.
//  Undocumented: should only be used inside Charm++.
/*forward*/ class Group;
class IrrGroup : public Chare {
  protected:
    CkGroupID thisgroup;
  public:
    IrrGroup(CkMigrateMessage *m): Chare(m) { }
    IrrGroup();
    virtual ~IrrGroup(); //<- needed for *any* child to have a virtual destructor

    virtual void pup(PUP::er &p);//<- pack/unpack routine
    virtual void ckJustMigrated(void);
    inline const CkGroupID &ckGetGroupID(void) const {return thisgroup;}
    inline CkGroupID CkGetGroupID(void) const {return thisgroup;}
    virtual int ckGetChareType() const;
    virtual char *ckDebugChareName();
    virtual int ckDebugChareID(char *, int);

    // Silly run-time type information
    virtual int isNodeGroup() { return 0; };
    virtual CmiBool isLocMgr(void){ return CmiFalse; }
    virtual CmiBool isArrMgr(void){ return CmiFalse; }
    virtual CmiBool isReductionMgr(void){ return CmiFalse; }
    static int isIrreducible(){ return 1;}
    virtual void flushStates() {}
		/*
			FAULT_EVAC
		*/
		virtual void evacuate(){};
		virtual void doneEvacuate(){};
    virtual void CkAddThreadListeners(CthThread tid, void *msg);
};

#define CBASE_PROXY_MEMBERS(CProxy_Derived) \
	typedef typename CProxy_Derived::local_t local_t; \
	typedef typename CProxy_Derived::index_t index_t; \
	typedef typename CProxy_Derived::proxy_t proxy_t; \
	typedef typename CProxy_Derived::element_t element_t; \
	CProxy_Derived thisProxy; 


/*Templated implementation of CBase_* classes.*/
template <class Parent,class CProxy_Derived>
class CBaseT : public Parent {
public:
	CBASE_PROXY_MEMBERS(CProxy_Derived)

	CBaseT(void) :Parent()  { thisProxy=this; }
	CBaseT(CkMigrateMessage *m) :Parent(m) { thisProxy=this; }
	void pup(PUP::er &p) {
		Parent::pup(p);
		p|thisProxy;
	}
};

/*Templated version of above for multiple (at least duplicate) inheritance:*/
template <class Parent1,class Parent2,class CProxy_Derived>
class CBaseT2 : public Parent1, public Parent2 {
public:
	CBASE_PROXY_MEMBERS(CProxy_Derived)

	CBaseT2(void) :Parent1(), Parent2()
		{ thisProxy = (Parent1 *)this; }
	CBaseT2(CkMigrateMessage *m) :Parent1(m), Parent2(m)
		{ thisProxy = (Parent1 *)this; } 
	void pup(PUP::er &p) {
		Parent1::pup(p);
		Parent2::pup(p);
		p|thisProxy;
	}

//These overloads are needed to prevent ambiguity for multiple inheritance:
	inline const CkChareID &ckGetChareID(void) const
		{return ((Parent1 *)this)->ckGetChareID();}
	static int isIrreducible(){ return (Parent1::isIrreducible() && Parent2::isIrreducible());}
	CHARM_INPLACE_NEW
};

/**************************** CkDelegateMgr **************************/

class CProxy;

/**
 Per-proxy data storage for delegation.  A CkDelegateMgr
 inherits from this class and adds his per-proxy data. 
 This class is reference counted.
*/
class CkDelegateData : public CkNoncopyable {
	int refcount; // reference count
public:
	CkDelegateData() :refcount(0) {}
	virtual ~CkDelegateData();
	
        //Child class constructor may have to set this.
        inline void reset() {
            refcount = 0;
        }
        
	/// Add a reference to this delegation data.  Just increments the refcount.
	///   Only CProxy should ever have to call this routine.
        /// Actually now the delegate manager calls it.
	inline void ref(void) {refcount++;}
	
	/// Remove our reference from this data.  If the refcount
	///  reaches 0, deletes the delegateData.
	///   Only CProxy should ever have to call this routine.
        /// Actually now the delegate manager calls it.
	inline void unref(void) {
		refcount--;
		if (refcount==0) delete this;
	}
};

/**
Message delegation support, where you send a message via
a proxy normally, but the message ends up routed via this
special delegateMgr group.

An "interface" class-- all delegated messages are routed via 
this class.  The default action is to deliver the message directly.
*/
class CkDelegateMgr : public IrrGroup {
  public:
    virtual ~CkDelegateMgr(); //<- so children can have virtual destructor
    virtual void ChareSend(CkDelegateData *pd,int ep,void *m,const CkChareID *c,int onPE);

    virtual void GroupSend(CkDelegateData *pd,int ep,void *m,int onPE,CkGroupID g);
    virtual void GroupBroadcast(CkDelegateData *pd,int ep,void *m,CkGroupID g);
    virtual void GroupSectionSend(CkDelegateData *pd,int ep,void *m,int nsid,CkSectionID *s);

    virtual void NodeGroupSend(CkDelegateData *pd,int ep,void *m,int onNode,CkNodeGroupID g);
    virtual void NodeGroupBroadcast(CkDelegateData *pd,int ep,void *m,CkNodeGroupID g);
    virtual void NodeGroupSectionSend(CkDelegateData *pd,int ep,void *m,int nsid,CkSectionID *s);

    virtual void ArrayCreate(CkDelegateData *pd,int ep,void *m,const CkArrayIndexMax &idx,int onPE,CkArrayID a);
    virtual void ArraySend(CkDelegateData *pd,int ep,void *m,const CkArrayIndexMax &idx,CkArrayID a);
    virtual void ArrayBroadcast(CkDelegateData *pd,int ep,void *m,CkArrayID a);
    virtual void ArraySectionSend(CkDelegateData *pd,int ep,void *m,int nsid,CkSectionID *s,int opts);
    virtual void initDelegateMgr(CProxy *proxy)  {}
    virtual CkDelegateData* ckCopyDelegateData(CkDelegateData *data) {
        data->ref();
        return data;
    } 
    
    /**
     Management of per-proxy data: pup this delegate's data.
     If p.isUnpacking, allocate and return a new set of delegate data.
     Never delete (or unref) the data-- the proxy will do that itself
        when it is required.
     The default implementation just ignores this call.
     
     A typical implementation that uses CkDelegateData might look like this:
     <code>
       myData *d=(myData *)pd;
       if (p.isUnpacking()) d=new myData();
       p|d->myField1;
       p|d->myField2;
       return d;
     </code>
    */
    virtual CkDelegateData *DelegatePointerPup(PUP::er &p,CkDelegateData *pd);
};


/**************************** Proxies **************************/

/*Message delegation support, where you send a message via
a proxy normally, but the message ends up routed via a 
special delegateMgr group.
*/
class CkDelegateMgr;

/** 
  A proxy is a local handle to a remote object.  This is the superclass
  of all proxies: CProxy_Array, CProxy_Group, etc. inherit from this class.
  
  Real proxies for user classes are generated by the .ci file translator charmxi
  and put in the generated .decl.h headers.
*/
class CProxy {
  private:
    CkDelegateMgr *delegatedMgr; // can be either a group or a nodegroup
    CkDelegateData *delegatedPtr; // private data for use by delegatedMgr.
  protected: //Never allocate CProxy's-- only subclass them.

#define CK_DELCTOR_PARAM CkDelegateMgr *dTo=NULL,CkDelegateData *dPtr=NULL
#define CK_DELCTOR_ARGS dTo,dPtr
#define CK_DELCTOR_CALL ckDelegatedTo(),ckDelegatedPtr()
/// Delegation constructor: used when building 
///   an element proxy from a collective proxy, like in "aProxy[i]".
    CProxy(CK_DELCTOR_PARAM)
	:delegatedMgr(dTo)
        {
            delegatedPtr = NULL;
            if(delegatedMgr != NULL && dPtr != NULL) 
                delegatedPtr = dTo->ckCopyDelegateData(dPtr);            
        }
  public:
    /// Copy constructor.  Only needed for delegated proxies.
    CProxy(const CProxy &src);
    /// Assignment operator.  Only needed for delegated proxies.
    CProxy& operator=(const CProxy &src);
    /// Destructor.  Only needed for delegated proxies. 
    ~CProxy() {
        if (delegatedPtr) delegatedPtr->unref();
    }
    
    /**
      Delegation allows a class, called a CkDelegateMgr, to 
      intercept calls made to this proxy for further processing.
      
      "ptr" is any delegator-specific data the CkDelegateMgr wants
      to associate with this proxy: the pointer is owned by this 
      proxy, but will be copied and pupped by calling delegator routines.
      
      This interface should only be used by library writers,
      not ordinary user code.
    */
    void ckDelegate(CkDelegateMgr *to,CkDelegateData *pd=NULL);
    
    /// Remove delegation from this proxy.
    void ckUndelegate(void);
    
    /// Return true if this proxy is delegated.
    int ckIsDelegated(void) const { return(delegatedMgr!=NULL);}
    
    /// Return the delegator of this proxy, to which the proxies' messages
    ///  are actually sent.
    inline CkDelegateMgr *ckDelegatedTo(void) const { return delegatedMgr; }
    
    /// Return the delegator's local data associated with this proxy.
    inline CkDelegateData *ckDelegatedPtr(void) const {return delegatedPtr;}
    
    /// Return the groupID of our delegator.
    ///   Note that this can be a GroupID or a NodeGroupID, so be careful!
    CkGroupID ckDelegatedIdx(void) const {
    	if (delegatedMgr) return delegatedMgr->CkGetGroupID();
	else {
	  CkGroupID gid; gid.setZero();
	  return gid;
	}
    }
    
    /// Pup the data for this proxy.  Only needed for delegated proxies.
    void pup(PUP::er &p);
};

PUPmarshall(CProxy)

/*These disambiguation macros are needed to support
  multiple inheritance in Chares (Groups, Arrays).
  They resolve ambiguous accessor calls to the parent "super".
  Because mutator routines need to change *all* the base
  classes, mutators are generated in xi-symbol.C.
*/
#define CK_DISAMBIG_CPROXY(super) \
    int ckIsDelegated(void) const {return super::ckIsDelegated();}\
    inline CkDelegateMgr *ckDelegatedTo(void) const {return super::ckDelegatedTo();}\
    inline CkDelegateData *ckDelegatedPtr(void) const {return super::ckDelegatedPtr();} \
    CkGroupID ckDelegatedIdx(void) const {return super::ckDelegatedIdx();}\



/*The base classes of each proxy type
*/
class CProxy_Chare : public CProxy {
  private:
    CkChareID _ck_cid;
  public:
    CProxy_Chare() {
#ifndef CMK_OPTIMIZE
	_ck_cid.onPE=0; _ck_cid.objPtr=0;
#endif
    }
#ifndef CMK_OPTIMIZE
    inline void ckCheck(void) const  {   //Make sure this proxy has a value
#if !CMK_FT_CHARE
	if (_ck_cid.objPtr==0)
		CkAbort("Error! This chare proxy has not been initialized!");
#endif
    }
#else
    inline void ckCheck() const {}
#endif
    CProxy_Chare(const CkChareID &c) : _ck_cid(c) {}
    CProxy_Chare(const Chare *c) : _ck_cid(c->ckGetChareID()) {}
    const CkChareID &ckGetChareID(void) const {return _ck_cid;}
    operator const CkChareID &(void) const {return ckGetChareID();}
    void ckSetChareID(const CkChareID &c) {_ck_cid=c;}
    void pup(PUP::er &p) {
    	CProxy::pup(p);
    	p(_ck_cid.onPE);
    	//Copy the pointer as straight bytes
    	p((char *)&_ck_cid.objPtr,sizeof(_ck_cid.objPtr));
    }
};
PUPmarshall(CProxy_Chare)

#define CK_DISAMBIG_CHARE(super) \
	CK_DISAMBIG_CPROXY(super) \
	inline void ckCheck(void) const {super::ckCheck();} \
	const CkChareID &ckGetChareID(void) const\
    	   {return super::ckGetChareID();} \
        operator const CkChareID &(void) const {return ckGetChareID();}

/******************* Reduction Declarations ****************/
//Silly: need the type of a reduction client here so it can be used by proxies.
//A clientFn is called on PE 0 when all contributions
// have been received and reduced.
//  param can be ignored, or used to pass any client-specific data you $
//  dataSize gives the size (in bytes) of the data array
//  data gives the reduced contributions--
//       it will be disposed of after this procedure returns.
typedef void (*CkReductionClientFn)(void *param,int dataSize,void *data);

/// Tiny utility class used by CkReductionClientAdaptor--
/// lets us keep backward compatability with the old C-style interface.
class CkReductionClientBundle : public CkCallback {
	CkReductionClientFn fn;
	void *param;
 public:
	static void callbackCfn(void *thisPtr,void *reductionMsg);
        CkReductionClientBundle(): fn(NULL), param(NULL) {}
	CkReductionClientBundle(CkReductionClientFn fn_,void *param_);
};
PUPbytes(CkReductionClientBundle)

#define CK_REDUCTION_CLIENT_DECL \
	void setReductionClient(CkReductionClientFn fn,void *param=NULL) const\
		{ ckSetReductionClient(fn,param); } \
	void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const \
		{ ckSetReductionClient(new CkReductionClientBundle(fn,param)); } \
	void ckSetReductionClient(CkCallback *cb) const;\

#define CK_REDUCTION_CLIENT_DEF(className,mgr) \
 	void className::ckSetReductionClient(CkCallback *cb) const \
		{ (mgr)->ckSetReductionClient(cb); }\

#define CK_REDUCTION_CLIENT_DISAMBIG(super) \
	inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const \
		{ super::setReductionClient(fn,param); } \
	inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const \
		{ super::ckSetReductionClient(fn,param); } \
	inline void ckSetReductionClient(CkCallback *cb) const \
		{ super::ckSetReductionClient(cb); }\

class CProxy_NodeGroup;
class CProxy_CkArrayReductionMgr;
class CProxy_Group : public CProxy {
  private:
    CkGroupID _ck_gid;

  public:
    CProxy_Group() {
#ifndef CMK_OPTIMIZE
	_ck_gid.setZero();
#endif
	//CkPrintf(" In CProxy_Group Constructor\n");
    }
    CProxy_Group(CkGroupID g,CK_DELCTOR_PARAM)
    	:CProxy(CK_DELCTOR_ARGS),_ck_gid(g) {
	//CkPrintf(" In CProxy_Group Constructor\n");
	}
    CProxy_Group(const IrrGroup *g)
        :CProxy(), _ck_gid(g->ckGetGroupID()) {
	//CkPrintf(" In CProxy_Group Constructor\n");
	}
/*    CProxy_Group(const NodeGroup *g)  //<- for compatability with NodeGroups
        :CProxy(), _ck_gid(g->ckGetGroupID()) {}*/

#ifndef CMK_OPTIMIZE
    inline void ckCheck(void) const {   //Make sure this proxy has a value
	if (_ck_gid.isZero())
		CkAbort("Error! This group proxy has not been initialized!");
    }
#else
    inline void ckCheck() const {}
#endif

    CkChareID ckGetChareID(void) const {
    	CkChareID ret;
    	ret.onPE=CkMyPe();
    	ret.objPtr=CkLocalBranch(_ck_gid);
    	return ret;
    }
    CkGroupID ckGetGroupID(void) const {return _ck_gid;}
    operator CkGroupID () const {return ckGetGroupID();}
    void ckSetGroupID(CkGroupID g) {_ck_gid=g;}
    void pup(PUP::er &p) {
    	CProxy::pup(p);
	p|_ck_gid;
    }
    CK_REDUCTION_CLIENT_DECL
};
PUPmarshall(CProxy_Group)
#define CK_DISAMBIG_GROUP(super) \
	CK_DISAMBIG_CPROXY(super) \
	inline void ckCheck(void) const {super::ckCheck();} \
	CkChareID ckGetChareID(void) const \
	   {return super::ckGetChareID();} \
	CkGroupID ckGetGroupID(void) const \
	   {return super::ckGetGroupID();} \
	operator CkGroupID () const { return ckGetGroupID(); } \
	CK_REDUCTION_CLIENT_DISAMBIG(super)\


class CProxyElement_Group : public CProxy_Group {
  private:
    int _onPE;
  public:
    CProxyElement_Group() { }
    CProxyElement_Group(CkGroupID g,int onPE,CK_DELCTOR_PARAM)
	: CProxy_Group(g,CK_DELCTOR_ARGS),_onPE(onPE) {}
    CProxyElement_Group(const IrrGroup *g)
        :CProxy_Group(g), _onPE(CkMyPe()) {}
    /*CProxyElement_Group(const NodeGroup *g)  //<- for compatability with NodeGroups
        :CProxy_Group(g), _onPE(CkMyPe()) {}*/

    int ckGetGroupPe(void) const {return _onPE;}
    void pup(PUP::er &p) {
    	CProxy_Group::pup(p);
    	p(_onPE);
    }
};
PUPmarshall(CProxyElement_Group)
#define CK_DISAMBIG_GROUP_ELEMENT(super) \
	CK_DISAMBIG_GROUP(super) \
	int ckGetGroupPe(void) const\
    	   {return super::ckGetGroupPe();} \


class CProxySection_Group : public CProxy_Group {
private:
  int _nsid;
  CkSectionID *_sid;
public:
  CProxySection_Group() { }
  CProxySection_Group(const CkGroupID &gid, const int *elems, const int nElems,CK_DELCTOR_PARAM)
      :CProxy_Group(gid,CK_DELCTOR_ARGS), _nsid(1) { _sid = new CkSectionID(gid, elems, nElems); }
  CProxySection_Group(const CProxySection_Group &cs,CK_DELCTOR_PARAM)
      :CProxy_Group(cs.ckGetGroupID(),CK_DELCTOR_ARGS), _nsid(cs._nsid) {
    if (_nsid == 1) _sid = new CkSectionID(cs.ckGetGroupID(), cs.ckGetElements(), cs.ckGetNumElements());
    else if (_nsid > 1) {
      _sid = new CkSectionID[_nsid];
      for (int i=0; i<_nsid; ++i) _sid[i] = cs._sid[i];
    } else _sid = NULL;
  }
  CProxySection_Group(const IrrGroup *g)
      :CProxy_Group(g), _nsid(0) {}
  CProxySection_Group(const int n, const CkGroupID *gid, const int **elems, const int *nElems,CK_DELCTOR_PARAM)
      :CProxy_Group(gid[0],CK_DELCTOR_ARGS), _nsid(n) {
    _sid = new CkSectionID[n];
    for (int i=0; i<n; ++i) _sid[i] = CkSectionID(gid[i], elems[i], nElems[i]);
  }
  
  ~CProxySection_Group() {
    if (_nsid == 1) delete _sid;
    else if (_nsid > 1) delete[] _sid;
  }
  
  CProxySection_Group &operator=(const CProxySection_Group &cs) {
    CProxy_Group::operator=(cs);
    _nsid = cs._nsid;
    if (_nsid == 1) _sid = new CkSectionID(*cs._sid);
    else if (_nsid > 1) {
      _sid = new CkSectionID[_nsid];
      for (int i=0; i<_nsid; ++i) _sid[i] = cs._sid[i];
    } else _sid = NULL;
    return *this;
  }
  
  //void ckSend(CkArrayMessage *m, int ep, int opts = 0) ;

  inline int ckGetNumSections() const {return _nsid;}
  inline CkSectionInfo &ckGetSectionInfo() {return _sid[0]._cookie;}
  inline CkSectionID *ckGetSectionID() {return _sid; }
  inline CkGroupID ckGetGroupIDn(int i) const {return (CkGroupID)_sid[i]._cookie.aid;}
  inline int *ckGetElements() const {return _sid[0].pelist;}
  inline int *ckGetElements(int i) const {return _sid[i].pelist;}
  inline int ckGetNumElements() const { return _sid[0].npes; }
  inline int ckGetNumElements(int i) const { return _sid[i].npes; }
  void pup(PUP::er &p) {
    CProxy_Group::pup(p);
    p | _nsid;
    if (p.isUnpacking()) {
      if (_nsid == 1) _sid = new CkSectionID;
      else if (_nsid > 1) _sid = new CkSectionID[_nsid];
      else _sid = NULL;
    }
    for (int i=0; i<_nsid; ++i) p | _sid[i];
  }
};
PUPmarshall(CProxySection_Group)
#define CK_DISAMBIG_GROUP_SECTION(super) \
    CK_DISAMBIG_GROUP(super) \
        inline int ckGetNumSections() const \
      { return super::ckGetNumSections(); } \
        inline CkSectionInfo &ckGetSectionInfo() \
      { return super::ckGetSectionInfo(); } \
        inline CkSectionID *ckGetSectionID() \
      { return super::ckGetSectionID(); } \
        inline CkGroupID ckGetGroupIDn(int i) const \
      { return super::ckGetGroupIDn(i); } \
        inline int *ckGetElements() const \
      { return super::ckGetElements(); } \
        inline int *ckGetElements(int i) const \
      { return super::ckGetElements(i); } \
        inline int ckGetNumElements() const \
      { return super::ckGetNumElements(); }  \
        inline int ckGetNumElements(int i) const \
      { return super::ckGetNumElements(i); }  \


/* These classes exist to provide chare indices for the basic
 chare types.*/
class CkIndex_Chare { public:
    static int __idx;//Fake chare index for registration
};
class CkIndex_ArrayBase { public:
    static int __idx;//Fake chare index for registration
};
class CkIndex_Group { public:
    static int __idx;//Fake chare index for registration
};

typedef CkIndex_Group CkIndex_NodeGroup;
typedef CkIndex_Group CkIndex_IrrGroup;


//typedef CProxy_Group CProxy_NodeGroup;
class CProxy_NodeGroup : public CProxy{

  private:
    CkGroupID _ck_gid;
  public:
    CProxy_NodeGroup() {
#ifndef CMK_OPTIMIZE
	_ck_gid.setZero();
#endif
	//CkPrintf("In CProxy_NodeGroup0 Constructor %d\n",CkLocalNodeBranch(_ck_gid));
    }
    CProxy_NodeGroup(CkGroupID g,CK_DELCTOR_PARAM)
    	:CProxy(CK_DELCTOR_ARGS),_ck_gid(g) {/*CkPrintf("In CProxy_NodeGroup2 Constructor %d\n",CkLocalNodeBranch(_ck_gid));*/}
    CProxy_NodeGroup(const IrrGroup *g)
        :CProxy(), _ck_gid(g->ckGetGroupID()) {/*CkPrintf("In CProxy_NodeGroup3 Constructor %d\n",CkLocalNodeBranch(_ck_gid));*/}
/*    CProxy_Group(const NodeGroup *g)  //<- for compatability with NodeGroups
        :CProxy(), _ck_gid(g->ckGetGroupID()) {}*/

#ifndef CMK_OPTIMIZE
    inline void ckCheck(void) const {   //Make sure this proxy has a value
	if (_ck_gid.isZero())
		CkAbort("Error! This group proxy has not been initialized!");
    }
#else
    inline void ckCheck() const {}
#endif

    CkChareID ckGetChareID(void) const {
    	CkChareID ret;
    	ret.onPE=CkMyPe();
    	ret.objPtr=CkLocalBranch(_ck_gid);
    	return ret;
    }
    CkGroupID ckGetGroupID(void) const {return _ck_gid;}
    operator CkGroupID () const {return ckGetGroupID();}
    void ckSetGroupID(CkGroupID g) {_ck_gid=g;}
    void pup(PUP::er &p) {
    	CProxy::pup(p);
	p | _ck_gid;
    }
    CK_REDUCTION_CLIENT_DECL

};

typedef CProxy_Group CProxy_IrrGroup;
typedef CProxyElement_Group CProxyElement_NodeGroup;
typedef CProxyElement_Group CProxyElement_IrrGroup;
typedef CProxySection_Group CProxySection_NodeGroup;
typedef CProxySection_Group CProxySection_IrrGroup;


//(CProxy_ArrayBase is defined in ckarray.h)

//Defines the actual "Group"
#include "ckreduction.h"

class CkQdMsg {
  public:
    void *operator new(size_t s) { return CkAllocMsg(0,(int)s,0); }
    void operator delete(void* ptr) { CkFreeMsg(ptr); }
    static void *alloc(int, size_t s, int*, int) {
      return CkAllocMsg(0,(int)s,0);
    }
    static void *pack(CkQdMsg *m) { return (void*) m; }
    static CkQdMsg *unpack(void *buf) { return (CkQdMsg*) buf; }
};

class CkThrCallArg {
  public:
    void *msg;
    void *obj;
    CkThrCallArg(void *m, void *o) : msg(m), obj(o) {}
};

extern void CkStartQD(const CkCallback& cb);
#define CkExitAfterQuiescence() CkStartQD(CkCallback(CkCallback::ckExit))


#if !CMK_MACHINE_PROGRESS_DEFINED
#define CkNetworkProgress() 
#define CkNetworkProgressAfter(p) 

#else
void CmiMachineProgressImpl();

#define CkNetworkProgress() {CpvAccess(networkProgressCount) ++; \
if(CpvAccess(networkProgressCount) >=  networkProgressPeriod)  \
    if (LBDatabaseObj()->getLBDB()->StatsOn() == 0) { \
        CmiMachineProgressImpl(); \
        CpvAccess(networkProgressCount) = 0; \
    } \
} \

#define CkNetworkProgressAfter(p) {CpvAccess(networkProgressCount) ++; \
if(CpvAccess(networkProgressCount) >=  p)  \
    if (LBDatabaseObj()->getLBDB()->StatsOn() == 0) { \
        CmiMachineProgressImpl(); \
        CpvAccess(networkProgressCount) = 0; \
    } \
} \

#endif


#ifdef _FAULT_MLOG_
#include "ckmessagelogging.h"
#endif
#include "ckmemcheckpoint.h"
#include "readonly.h"
#include "ckarray.h"
#include "ckstream.h"
#include "CkFutures.decl.h"
#include "charisma.h"
#include "tempo.h"
#include "waitqd.h"
#include "sdag.h"
#include "ckcheckpoint.h"
#include "ckevacuation.h"
#include "ckarrayreductionmgr.h"
#include "trace.h"
#include "envelope.h"






CkMarshallMsg *CkAllocateMarshallMsgNoninline(int size,const CkEntryOptions *opts);
inline CkMarshallMsg *CkAllocateMarshallMsg(int size,const CkEntryOptions *opts=NULL)
{
	if (opts==NULL) {
	  CkMarshallMsg *newMemory = new (size,0)CkMarshallMsg;
	  setMemoryTypeMessage(UsrToEnv(newMemory));
	  return newMemory;
	}
	else return CkAllocateMarshallMsgNoninline(size,opts);
}







template <typename T> 
inline T *CkAllocateMarshallMsgT(int size,const CkEntryOptions *opts) 
{ 
  int priobits = 0; 
  if (opts!=NULL) priobits = opts->getPriorityBits(); 
  //Allocate the message 
  T *m=new (size,priobits)T; 
  //Copy the user's priority data into the message 
  envelope *env=UsrToEnv(m); 
  setMemoryTypeMessage(env); 
  if (opts!=NULL) { 
    CmiMemcpy(env->getPrioPtr(),opts->getPriorityPtr(),env->getPrioBytes()); 
    //Set the message's queueing type 
    env->setQueueing((unsigned char)opts->getQueueing()); 
  } 
  return m; 
} 





/************************** Debugging Utilities **************/

//For debugging: convert given index to a string (NOT threadsafe)
static const char *idx2str(const CkArrayIndex &ind) {
  static char retBuf[80];
  retBuf[0]=0;
  if (ind.dimension <= 3) {
    for (int i=0;i<ind.nInts;i++) {
      if (i>0) strcat(retBuf,";");
      sprintf(&retBuf[strlen(retBuf)],"%d",ind.data()[i]);
    }
  } else {
    const short int *idx = (const short int*)ind.data();
    for (int i=0;i<ind.dimension;i++) {
      if (i>0) strcat(retBuf,";");
      sprintf(&retBuf[strlen(retBuf)],"%hd",idx[i]);
    }
  }
  return retBuf;
}

static const char *idx2str(const ArrayElement *el) {
  return idx2str(el->thisIndexMax);
}



class CkConditional {
  int refcount;
public:
  CkConditional() : refcount(1) { }
  virtual ~CkConditional() { }
    /*
  void* operator new(size_t s) {
    return malloc(s);
  }
    */
  void deallocate() {
    CkPrintf("CkConditional::delete %d\n",refcount);
    if (--refcount == 0) {
      //((CkConditional*)p)->~CkConditional();
      delete this;
    }
  }
  void copyreference(void) {
    ++refcount;
  }
};


#endif



