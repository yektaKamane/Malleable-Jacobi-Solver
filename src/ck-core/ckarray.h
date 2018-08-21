/* Generalized Chare Arrays

These classes implement Chare Arrays.
These are dynamic (i.e. allowing insertion
and deletion) collections of ordinary Chares
indexed by arbitrary runs of bytes.

The general structure is:

CkArray is the "array manager" Group, or BOC--
it creates, keeps track of, and cares for all the
array elements on this PE (i.e.. "local" elements).
It does so using a hashtable of CkArrayRec objects--
there's an entry for each local, home-here, and recently
communicated remote array elements.

CkArrayElement is the type of the array
elements (a subclass of Chare).

CkArrayIndex is an arbitrary run of bytes,
used to index into the CkArray hashtable.

Converted from 1-D arrays 2/27/2000 by
Orion Sky Lawlor, olawlor@acm.org
*/
#ifndef __CKARRAY_H
#define __CKARRAY_H

#include "cklocation.h"
#include "ckmulticast.h"
#include "ckmemcheckpoint.h" // for CkArrayCheckPTReqMessage
#include "ckarrayindex.h"

/***********************************************************
	Utility defines, includes, etc.
*/
extern void _registerCkArray(void);
CpvExtern (int ,serializer);

#if (defined(_FAULT_MLOG_) || defined(_FAULT_CAUSAL_))
#define _MLOG_BCAST_TREE_ 1
#define _MLOG_BCAST_BFACTOR_ 8
#endif

/** This flag is true when in the system there is anytime migration, false when
 *  the user code guarantees that no migration happens except during load balancing
 *  (in which case it can only happen between AtSync and ResumeFromSync). */
extern bool _isAnytimeMigration;

/**
  Array elements are only inserted at construction
 */
extern bool _isStaticInsertion;

/** This flag is true when users are sure there is at least one charm array element
 *  per processor. In such case, when doing reduction on the array, the children
 *  don't need to be notified that reduction starts
 */
extern bool _isNotifyChildInRed;

/**
\addtogroup CkArray
\brief Migratable Chare Arrays: user-visible classes.

All these classes are defined in ckarray.C.
*/
/*@{*/
/*@}*/

//This simple arrayListener just prints each event to stdout:
class CkVerboseListener : public CkArrayListener {
 public:
  CkVerboseListener(void);
  CkVerboseListener(CkMigrateMessage *m):CkArrayListener(m) {}
  PUPable_decl(CkVerboseListener);

  virtual void ckRegister(CkArray *arrMgr,int dataOffset_);
  virtual void ckBeginInserting(void);
  virtual void ckEndInserting(void);

  virtual void ckElementStamp(int *eltInfo);
  virtual void ckElementCreating(ArrayElement *elt);
  virtual bool ckElementCreated(ArrayElement *elt);
  virtual void ckElementDied(ArrayElement *elt);

  virtual void ckElementLeaving(ArrayElement *elt);
  virtual bool ckElementArriving(ArrayElement *elt);
};

/**
\addtogroup CkArray
*/
/*@{*/

/*********************** Proxy Support ************************/
//Needed by CBase_ArrayElement
class ArrayBase { /*empty*/ };
/*forward*/ class ArrayElement;

/**
 * This class is a wrapper around a CkArrayIndex and ArrayID,
 *  used by array element proxies.  This makes the translator's
 *  job simpler, and the translated code smaller.
 */
class CProxy_ArrayBase :public CProxy {
private:
	CkArrayID _aid;
public:
	CProxy_ArrayBase() {
#if CMK_ERROR_CHECKING
		_aid.setZero();
#endif
	}
	CProxy_ArrayBase(const CkArrayID &aid,CK_DELCTOR_PARAM)
		:CProxy(CK_DELCTOR_ARGS), _aid(aid) { }
        CProxy_ArrayBase(const CkArrayID &aid)
                :CProxy(), _aid(aid) { }
	CProxy_ArrayBase(const ArrayElement *e);
	CProxy_ArrayBase(const CProxy_ArrayBase &cs): CProxy(cs), _aid(cs.ckGetArrayID()) {}

	bool operator==(const CProxy_ArrayBase& other) {
		return ckGetArrayID() == other.ckGetArrayID();
	}

#if CMK_ERROR_CHECKING
	inline void ckCheck(void) const{  //Make sure this proxy has a value
	  if (_aid.isZero())
		CkAbort("Error! This array proxy has not been initialized!");
        }
#else
	inline void ckCheck(void) const {}
#endif

	static CkArrayID ckCreateEmptyArray(CkArrayOptions opts);
        static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts);
	static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts);

	void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx);
	void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const;
	CkArrayID ckGetArrayID(void) const { return _aid; }
	CkArray *ckLocalBranch(void) const { return _aid.ckLocalBranch(); }
	CkLocMgr *ckLocMgr(void) const;
	inline operator CkArrayID () const {return ckGetArrayID();}
	unsigned int numLocalElements() const { return ckLocMgr()->numLocalElements(); }

	void doneInserting(void);
	void beginInserting(void);

	CK_REDUCTION_CLIENT_DECL

	void pup(PUP::er &p);
};
PUPmarshall(CProxy_ArrayBase)

class CProxyElement_ArrayBase:public CProxy_ArrayBase {
private:
	CkArrayIndex _idx;//<- our element's array index
public:
	CProxyElement_ArrayBase() { }
	CProxyElement_ArrayBase(const CkArrayID &aid,
		const CkArrayIndex &idx,CK_DELCTOR_PARAM)
		:CProxy_ArrayBase(aid,CK_DELCTOR_ARGS), _idx(idx) { }
        CProxyElement_ArrayBase(const CkArrayID &aid, const CkArrayIndex &idx)
                :CProxy_ArrayBase(aid), _idx(idx) { }
	CProxyElement_ArrayBase(const ArrayElement *e);

	bool operator==(const CProxyElement_ArrayBase& other) {
		return ckGetArrayID() == other.ckGetArrayID() &&
				ckGetIndex() == other.ckGetIndex();
	}

	void ckInsert(CkArrayMessage *m,int ctor,int onPe);
	void ckSend(CkArrayMessage *m, int ep, int opts = 0) const;
//      static void ckSendWrapper(void *me, void *m, int ep, int opts = 0);
      static void ckSendWrapper(CkArrayID _aid, CkArrayIndex _idx, void *m, int ep, int opts);
	void *ckSendSync(CkArrayMessage *m, int ep) const;
	const CkArrayIndex &ckGetIndex() const {return _idx;}

	ArrayElement *ckLocal(void) const;
	void pup(PUP::er &p);
};
PUPmarshall(CProxyElement_ArrayBase)


#define _AUTO_DELEGATE_MCASTMGR_ON_ 1

class CProxySection_ArrayBase:public CProxy_ArrayBase {
private:
  std::vector<CkSectionID> _sid;
public:
	CProxySection_ArrayBase() =default;
	CProxySection_ArrayBase(const CkArrayID &aid,
		const CkArrayIndex *elems, const int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR);
	CProxySection_ArrayBase(const CkArrayID &aid,
		const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR);
	CProxySection_ArrayBase(const CkArrayID &aid,
		const CkArrayIndex *elems, const int nElems, CK_DELCTOR_PARAM)
		:CProxy_ArrayBase(aid,CK_DELCTOR_ARGS) { _sid.emplace_back(aid, elems, nElems); }
	CProxySection_ArrayBase(const CkArrayID &aid,
		const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM)
		:CProxy_ArrayBase(aid,CK_DELCTOR_ARGS) { _sid.emplace_back(aid, elems); }
	CProxySection_ArrayBase(const CkSectionID &sid)
		:CProxy_ArrayBase(sid._cookie.get_aid()) { _sid.emplace_back(sid); }
	CProxySection_ArrayBase(const CkSectionID &sid, CK_DELCTOR_PARAM)
		:CProxy_ArrayBase(sid._cookie.get_aid(), CK_DELCTOR_ARGS) { _sid.emplace_back(sid); }
        CProxySection_ArrayBase(const CProxySection_ArrayBase &cs)
		:CProxy_ArrayBase(cs) {
        _sid.resize(cs._sid.size());
        for (int i=0; i<_sid.size(); ++i) {
          _sid[i] = cs._sid[i];
        }
    }
        CProxySection_ArrayBase(const CProxySection_ArrayBase &cs, CK_DELCTOR_PARAM)
		:CProxy_ArrayBase(cs.ckGetArrayID(),CK_DELCTOR_ARGS) {
        _sid.resize(cs._sid.size());
        for (int i=0; i<_sid.size(); ++i) {
          _sid[i] = cs._sid[i];
        }
    }
    CProxySection_ArrayBase(const int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor=USE_DEFAULT_BRANCH_FACTOR);
    CProxySection_ArrayBase(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor=USE_DEFAULT_BRANCH_FACTOR);
    CProxySection_ArrayBase(const int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems,CK_DELCTOR_PARAM)
        :CProxy_ArrayBase(aid[0],CK_DELCTOR_ARGS) {
        _sid.resize(n);
        for (int i=0; i<_sid.size(); ++i) {
          _sid[i] = CkSectionID(aid[i], elems[i], nElems[i]);
        }
    }
    CProxySection_ArrayBase(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM)
        :CProxy_ArrayBase(aid[0],CK_DELCTOR_ARGS) {
        _sid.resize(aid.size());
        for (int i=0; i<_sid.size(); ++i) {
          _sid[i] = CkSectionID(aid[i], elems[i]);
        }
    }

    ~CProxySection_ArrayBase() =default;

    CProxySection_ArrayBase &operator=(const CProxySection_ArrayBase &cs) {
      CProxy_ArrayBase::operator=(cs);
      _sid.resize(cs._sid.size());
      for (int i=0; i<_sid.size(); ++i) {
        _sid[i] = cs._sid[i];
      }
      return *this;
    }
    
	void ckAutoDelegate(int opts=1);
 	using CProxy_ArrayBase::setReductionClient; //compilation error o/w
	void setReductionClient(CkCallback *cb); 
	void resetSection();

	void ckSectionDelegate(CkDelegateMgr *d, int opts=1) {
           ckDelegate(d);
	   if(opts==1)
              d->initDelegateMgr(this);
        }
//	void ckInsert(CkArrayMessage *m,int ctor,int onPe);
	void ckSend(CkArrayMessage *m, int ep, int opts = 0) ;

//	ArrayElement *ckLocal(void) const;
    inline int ckGetNumSubSections() const { return _sid.size(); }
	inline CkSectionInfo &ckGetSectionInfo() {return _sid[0]._cookie;}
	inline CkSectionID *ckGetSectionIDs() {return _sid.data();}
	inline CkSectionID &ckGetSectionID() {return _sid[0];}
	inline CkSectionID &ckGetSectionID(int i) {return _sid[i];}
	inline CkArrayID ckGetArrayIDn(int i) const {return _sid[i]._cookie.get_aid();}
    inline CkArrayIndex *ckGetArrayElements() const {return const_cast<CkArrayIndex *>(_sid[0]._elems.data());}
    inline CkArrayIndex *ckGetArrayElements(int i) const {return const_cast<CkArrayIndex *>(_sid[i]._elems.data());}
    inline int ckGetNumElements() const { return _sid[0]._elems.size(); }
	inline int ckGetNumElements(int i) const { return _sid[i]._elems.size(); }
	inline int ckGetBfactor() const { return _sid[0].bfactor; }
	void pup(PUP::er &p);
};
PUPmarshall(CProxySection_ArrayBase)

//Simple C-like API:
void CkSetMsgArrayIfNotThere(void *msg);
void CkSendMsgArray(int entryIndex, void *msg, CkArrayID aID, const CkArrayIndex &idx, int opts=0);
void CkSendMsgArrayInline(int entryIndex, void *msg, CkArrayID aID, const CkArrayIndex &idx, int opts=0);
void CkBroadcastMsgArray(int entryIndex, void *msg, CkArrayID aID, int opts=0);
void CkBroadcastMsgSection(int entryIndex, void *msg, CkSectionID sID, int opts=    0);
/************************ Array Element *********************/
/**
 *An array element is a chare that lives inside the array.
 *Unlike regular chares, array elements can migrate from one
 *Pe to another.  Each element has a unique index.
 */

class ArrayElement : public CkMigratable
{
  friend class CkArray;

  friend class CkArrayListener;
  int numInitialElements; // Number of elements created by ckNew(numElements)
  void initBasics(void);
#ifdef _PIPELINED_ALLREDUCE_
AllreduceMgr * allredMgr; // for allreduce
#endif
public:
  ArrayElement(void);
  ArrayElement(CkMigrateMessage *m);
  virtual ~ArrayElement();

/// Pack/unpack routine (called before and after migration)
  void pup(PUP::er &p);

//Overridden functions:
  /// Called by the system just before and after migration to another processor:
  virtual void ckAboutToMigrate(void);
  virtual void ckJustMigrated(void);
  
  virtual void ckJustRestored(void);
  
  virtual void ckDestroy(void);
  virtual char *ckDebugChareName(void);
  virtual int ckDebugChareID(char*, int);

  /// Synonym for ckMigrate
  inline void migrateMe(int toPe) {ckMigrate(toPe);}

#ifdef _PIPELINED_ALLREDUCE_
	void contribute2(CkArrayIndex myIndex, int dataSize,const void *data,CkReduction::reducerType type,
			   const CkCallback &cb,CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1);
	void contribute2(int dataSize,const void *data,CkReduction::reducerType type, 
					CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1); 
	void contribute2(int dataSize,const void *data,CkReduction::reducerType type, 
					const CkCallback &cb,CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1); 
	void contribute2(CkReductionMsg *msg); 
	void contribute2(const CkCallback &cb,CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1);
	void contribute2(CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1);
#else
  CK_REDUCTION_CONTRIBUTE_METHODS_DECL
#endif
	// for _PIPELINED_ALLREDUCE_, assembler entry method
	inline void defrag(CkReductionMsg* msg);
  inline const CkArrayID &ckGetArrayID(void) const {return thisArrayID;}
  inline ck::ObjID ckGetID(void) const { return ck::ObjID(thisArrayID, myRec->getID()); }

  inline int ckGetArraySize(void) const { return numInitialElements; }
protected:
  CkArray *thisArray;//My source array
  CkArrayID thisArrayID;//My source array's ID

  //More verbose form of abort
  virtual void CkAbort(const char *str) const;

private:
//Array implementation methods:
  int listenerData[CK_ARRAYLISTENER_MAXLEN];

#if CMK_MEM_CHECKPOINT
friend class CkMemCheckPT;
friend class CkLocMgr;
protected:
  int budPEs[2];
private:
  void init_checkpt();
#endif
public:
	void inmem_checkpoint(CkArrayCheckPTReqMessage *m);
	void recvBroadcast(CkMessage *);

#if CMK_GRID_QUEUE_AVAILABLE
public:
  int grid_queue_interval;
  int grid_queue_threshold;
  int msg_count;
  int msg_count_grid;
  int border_flag;
#endif
};
inline int *CkArrayListener::ckGetData(ArrayElement *el) const
  {return &el->listenerData[dataOffset];}

/**An ArrayElementT is a utility class where you are
 * constrained to a "thisIndex" of some fixed-sized type T.
 */
template <class T>
class ArrayElementT : public ArrayElement
{
public:
  using array_index_t = T;

  ArrayElementT(void): thisIndex(*(const T *)thisIndexMax.data()) {
#if (defined(_FAULT_MLOG_) || defined(_FAULT_CAUSAL_))     
        mlogData->objID.data.array.idx=thisIndexMax;
#endif
}
#ifdef _PIPELINED_ALLREDUCE_
	void contribute(int dataSize,const void *data,CkReduction::reducerType type,
						CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1)
	{
		contribute2( dataSize,data, type, userFlag);

	}
	void contribute(int dataSize,const void *data,CkReduction::reducerType type,
						const CkCallback &cb,CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1)
	{
		contribute2((CkArrayIndex)(thisIndex) ,dataSize, data, type,   cb, userFlag);
	}
	void contribute(CkReductionMsg *msg) 
	{
		contribute2(msg);
	}
	void contribute(const CkCallback &cb,CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1)
	{
		contribute2(cb ,userFlag);
	}
	void contribute(CMK_REFNUM_TYPE userFlag=(CMK_REFNUM_TYPE)-1)
	{
		contribute2(userFlag);
	}
#endif
  ArrayElementT(CkMigrateMessage *msg)
	:ArrayElement(msg),
	thisIndex(*(const T *)thisIndexMax.data()) {}

  const T thisIndex;/// Object array index
};

typedef ArrayElementT<CkIndex1D> ArrayElement1D;
typedef ArrayElementT<CkIndex2D> ArrayElement2D;
typedef ArrayElementT<CkIndex3D> ArrayElement3D;
typedef ArrayElementT<CkIndex4D> ArrayElement4D;
typedef ArrayElementT<CkIndex5D> ArrayElement5D;
typedef ArrayElementT<CkIndex6D> ArrayElement6D;
typedef ArrayElementT<CkIndexMax> ArrayElementMax;

extern void (*ArrayMsgRecvExtCallback)(int, int, int *, int, int, char *, int);
extern int (*ArrayElemLeaveExt)(int, int, int *, char**, int);
extern void (*ArrayElemJoinExt)(int, int, int *, int, char*, int);
extern void (*ArrayResumeFromSyncExtCallback)(int, int, int *);

class ArrayElemExt: public ArrayElement {
private:
  int ctorEpIdx;

public:
  ArrayElemExt(void *impl_msg);
  ArrayElemExt(CkMigrateMessage *m) { delete m; }

  static void __ArrayElemExt(void *impl_msg, void *impl_obj_void) {
    new (impl_obj_void) ArrayElemExt(impl_msg);
  }

  static void __entryMethod(void *impl_msg, void *impl_obj_void) {
    //fprintf(stderr, "ArrayElemExt:: Entry method invoked\n");
    ArrayElemExt *e = static_cast<ArrayElemExt *>(impl_obj_void);
    CkMarshallMsg *impl_msg_typed = (CkMarshallMsg *)impl_msg;
    char *impl_buf = impl_msg_typed->msgBuf;
    PUP::fromMem implP(impl_buf);
    int msgSize; implP|msgSize;
    int ep; implP|ep;
    int dcopy_start; implP|dcopy_start;
    ArrayMsgRecvExtCallback(((CkGroupID)e->thisArrayID).idx,
                            int(e->thisIndexMax.getDimension()),
                            e->thisIndexMax.data(), ep, msgSize, impl_buf+(3*sizeof(int)),
                            dcopy_start);
  }

  static void __AtSyncEntryMethod(void *impl_msg, void *impl_obj_void) {
    ArrayElemExt *e = static_cast<ArrayElemExt *>(impl_obj_void);
    //printf("ArrayElementExt:: calling AtSync elem->usesAtSync=%d\n", e->usesAtSync);
    e->AtSync();
    if (UsrToEnv(impl_msg)->isVarSysMsg() == 0) CkFreeSysMsg(impl_msg);
  }

  static void __CkMigrateMessage(void *impl_msg, void *impl_obj_void) {
    //printf("ArrayElemExt:: Migration constructor invoked\n");
    call_migration_constructor<ArrayElemExt> c = impl_obj_void;
    c((CkMigrateMessage*)impl_msg);
  }

  // NOTE this assumes that calls to pup are due to array element migration
  // not sure if it's always going to be the case
  void pup(PUP::er &p) {
    // because this is not generated code, looks like I need to make explicit call
    // to parent pup method, otherwise fields in parent class like usesAtSync will
    // not be pupped!
    ArrayElement::pup(p);
    int nDims = thisIndexMax.getDimension();
    int aid = ((CkGroupID)thisArrayID).idx;
    int data_size;
    if (!p.isUnpacking()) { // packing or sizing
      char *msg;
      data_size = ArrayElemLeaveExt(aid, nDims, thisIndexMax.data(), &msg, p.isSizing());
      p | data_size;
      p | ctorEpIdx;
      p(msg, data_size);
    } else {
      p | data_size;
      p | ctorEpIdx;
      PUP::fromMem *p_mem = (PUP::fromMem *)&p;
      ArrayElemJoinExt(aid, nDims, thisIndexMax.data(), ctorEpIdx, p_mem->get_current_pointer(), data_size);
      p_mem->advance(data_size);
    }
  }

  void ResumeFromSync() {
    ArrayResumeFromSyncExtCallback(((CkGroupID)thisArrayID).idx,
                            int(thisIndexMax.getDimension()),
                            thisIndexMax.data());
  }
};

/*@}*/


/*********************** Array Manager BOC *******************/
/**
\addtogroup CkArrayImpl
*/
/*@{*/

#include "CkArray.decl.h"

void CkSendAsyncCreateArray(int ctor, CkCallback cb, CkArrayOptions opts, void *ctorMsg);

class CkArrayCreatedMsg : public CMessage_CkArrayCreatedMsg {
public:
  CkArrayID aid;
  CkArrayCreatedMsg(CkArrayID _aid) : aid(_aid) { }
};

class CkArrayBroadcaster;
class CkArrayReducer;

void _ckArrayInit(void);

class CkArray : public CkReductionMgr {
  friend class ArrayElement;
  friend class CProxy_ArrayBase;
  friend class CProxyElement_ArrayBase;

  CkMagicNumber<ArrayElement> magic; //To detect heap corruption
  CkLocMgr *locMgr;
  CkGroupID locMgrID;
  CkGroupID mCastMgrID;
  bool sectionAutoDelegate;
  CkCallback initCallback;
  CProxy_CkArray thisProxy;
  // Separate mapping and storing the element pointers to speed iteration in broadcast
  std::map<CmiUInt8, unsigned int> localElems;
  std::vector<CkMigratable *> localElemVec;
#if (defined(_FAULT_MLOG_) || defined(_FAULT_CAUSAL_))
    int *children;
    int numChildren;
#endif
private:
  bool stableLocations;

public:
//Array Creation:
  CkArray(CkArrayOptions &&c, CkMarshalledMessage &&initMsg);
  CkArray(CkMigrateMessage *m);
  ~CkArray();
  CkGroupID &getGroupID(void) {return thisgroup;}
  CkGroupID &getmCastMgr(void) {return mCastMgrID;}
  bool isSectionAutoDelegated(void) {return sectionAutoDelegate;}

//Access & information routines
  inline CkLocMgr *getLocMgr(void) {return locMgr;}
  inline const CkArrayIndex &getNumInitial(void) const {return numInitial;}
  inline int homePe(const CkArrayIndex &idx) const {return locMgr->homePe(idx);}
  inline int procNum(const CkArrayIndex &idx) const {return locMgr->procNum(idx);}

  /// Return the last known processor for this array index.
  /// Valid for any possible array index.
  inline int lastKnown(const CkArrayIndex &idx) const
	  {return locMgr->lastKnown(idx);}
  /// Deliver message to this element (directly if local)
  /// doFree if is local
  inline void deliver(CkMessage *m, const CkArrayIndex &idx, CkDeliver_t type,int opts=0)
  { locMgr->sendMsg((CkArrayMessage*)m, thisgroup, idx, type, opts); }
  inline int deliver(CkArrayMessage *m, CkDeliver_t type)
  { return locMgr->deliverMsg(m, thisgroup, m->array_element_id(), NULL, type); }
  /// Fetch a local element via its ID (return NULL if not local)
  inline ArrayElement *lookup(const CmiUInt8 id) { return (ArrayElement*) getEltFromArrMgr(id); }
  /// Fetch a local element via its index (return NULL if not local)
  inline ArrayElement *lookup(const CkArrayIndex &idx) { 
    CmiUInt8 id;
    if (locMgr->lookupID(idx,id)) {
      return (ArrayElement*) getEltFromArrMgr(id);
    } else {
      return NULL;
    }
  }

  virtual CkMigratable* getEltFromArrMgr(const CmiUInt8 id) {
    std::map<CmiUInt8, unsigned int>::iterator itr = localElems.find(id);
    return ( itr == localElems.end() ? NULL : localElemVec[itr->second] );
  }
  virtual void putEltInArrMgr(const CmiUInt8 id, CkMigratable* elt)
  {
    localElems[id] = localElemVec.size();
    localElemVec.push_back(elt);
  }
  virtual void eraseEltFromArrMgr(const CmiUInt8 id)
  { localElems.erase(id); }

  void deleteElt(const CmiUInt8 id) {
    std::map<CmiUInt8, unsigned int>::iterator itr = localElems.find(id);
    if (itr != localElems.end()) {
      unsigned int offset = itr->second;
      localElems.erase(itr);
      delete localElemVec[offset];

      if (offset != localElemVec.size() - 1) {
        CkMigratable *moved = localElemVec[localElemVec.size()-1];
        localElemVec[offset] = moved;
        localElems[moved->ckGetID()] = offset;
      }

      localElemVec.pop_back();
    }
  }

//Creation:
  /// Create-after-migrate:
  /// Create an uninitialized element after migration
  ///  The element's constructor will be called immediately after.
  virtual CkMigratable *allocateMigrated(int elChareType, CkElementCreation_t type);
  void stampListenerData(CkMigratable *elt);

  /// Prepare creation message:
  void prepareCtorMsg(CkMessage *m, int listenerData[CK_ARRAYLISTENER_MAXLEN]);

  int findInitialHostPe(const CkArrayIndex &idx, int proposedPe);

  /// Create initial array elements:
  virtual void insertInitial(const CkArrayIndex &idx,void *ctorMsg);
  virtual void doneInserting(void);
  virtual void beginInserting(void);
  void remoteDoneInserting(void);
  void remoteBeginInserting(void);
  void initDone(void);

  /// Create manually:
  bool insertElement(CkArrayMessage *, const CkArrayIndex &idx, int listenerData[CK_ARRAYLISTENER_MAXLEN]);
  void insertElement(CkMarshalledMessage &&, const CkArrayIndex &idx, int listenerData[CK_ARRAYLISTENER_MAXLEN]);

/// Demand-creation:
  /// Demand-create an element at this index on this processor
  void demandCreateElement(const CkArrayIndex &idx, int ctor, CkDeliver_t type);

/// Broadcast communication:
  void sendBroadcast(CkMessage *msg);
  void recvBroadcast(CkMessage *msg);
  void sendExpeditedBroadcast(CkMessage *msg);
  void recvExpeditedBroadcast(CkMessage *msg) { recvBroadcast(msg); }
  void recvBroadcastViaTree(CkMessage *msg);

  /// Whole array destruction, including all elements and the group itself
  void ckDestroy();

  void pup(PUP::er &p);
  void ckJustMigrated(void){ doneInserting(); }

  virtual bool isArrMgr(void) {return true;}

private:
  CkArrayIndex numInitial;/// Number of initial array elements
  bool isInserting;/// Are we currently inserting elements?
  int numPesInited;

/// Allocate space for a new array element
  ArrayElement *allocate(int elChareType, CkMessage *msg, bool fromMigration, int *listenerData);

//Spring cleaning
  void springCleaning(void);
  static void staticSpringCleaning(void *forWhom,double curWallTime);
  void setupSpringCleaning();
  int springCleaningCcd;

//ArrayListeners:
  CkPupAblePtrVec<CkArrayListener> listeners;
  int listenerDataOffset;
 public:
  void addListener(CkArrayListener *l) {
    l->ckRegister(this,listenerDataOffset);
    listenerDataOffset+=l->ckGetLen();
    listeners.push_back(l);
    if (listenerDataOffset>CK_ARRAYLISTENER_MAXLEN)
      CkAbort("Too much array listener data!\n"
	      "You'll have to either use fewer array listeners, or increase the compile-time\n"
	      "constant CK_ARRAYLISTENER_MAXLEN!\n");
  }
 private:

  CkArrayReducer *reducer; //Read-only copy of default reducer
  CkArrayBroadcaster *broadcaster; //Read-only copy of default broadcaster
public:
  void flushStates();

#if (defined(_FAULT_MLOG_) || defined(_FAULT_CAUSAL_))
	// the mlogft only support 1D arrays, then returning the number of elements in the first dimension
	virtual int numberReductionMessages(){CkAssert(CkMyPe() == 0);return numInitial.data()[0];}
	void broadcastHomeElements(void *data,CkLocRec *rec,CkArrayIndex *index);
	static void staticBroadcastHomeElements(CkArray *arr,void *data,CkLocRec *rec,CkArrayIndex *index);
#endif

        static bool isIrreducible() { return true; }
};

// Maintain old name of former parent class for backwards API compatibility
// with usage in maps' populateInitial()
typedef CkArray CkArrMgr;

/*@}*/

#endif
