/*
Pack/UnPack Library for UIUC Parallel Programming Lab
Orion Sky Lawlor, olawlor@uiuc.edu, 4/5/2000

This library allows you to easily pack an array, structure,
or object into a memory buffer or disk file, and then read
the object back later.  The library can also handle translating
between different machine representations for integers and floats.

Typically, the user has to write separate functions for buffer 
sizing, pack to memory, unpack from memory, pack to disk, and 
unpack from disk.  These functions all perform the exact same function--
namely, they list the members of the array, struct, or object.
Further, all the functions must agree, or the unpacked data will 
be garbage.  This library allows the user to write *one* function,
pup, which will perform all needed packing/unpacking.

A simple example is:
class foo {
 private:
  CmiBool isBar;
  int x;
  char y;
  unsigned long z;
  float q[3];
 public:
  ...
  void pup(PUP::er &p) {
    p(isBar);
    p(x);p(y);p(z);
    p(q,3);
  }
};

A more complex example is:
class bar {
 private:
  foo f;
  int nArr;//Length of array below
  double *arr;//Heap-allocated array
 public:
  ...
  
  void pup(PUP::er &p) {
    f.pup(p);
    p(nArr);
    if (p.isUnpacking())
      arr=new double[nArr];
    p(arr,nArr);
  }
};
*/

#ifndef __CK_PUP_H
#define __CK_PUP_H

#include <stdio.h> /*<- for "FILE *" */

#ifndef __cplusplus
#error "Use pup_c.h for C programs-- pup.h is for C++ programs"
#endif

#ifndef CHARM_H
#  include <converse.h> // <- for CmiBool, CMK_* defines
#endif

//We need CkMigrateMessage only to distinguish the migration
// constructor from all other constructors-- the type
// itself has no meaningful fields.
typedef struct {int is_only_a_name;} CkMigrateMessage;

namespace PUP {

#if CMK_LONG_LONG_DEFINED
#define CMK_PUP_LONG_LONG long long
#elif CMK___int64_DEFINED
#define CMK_PUP_LONG_LONG __int64
#endif

 
//Item data types-- these are used to do byte swapping, etc.
typedef enum {
//(this list must exactly match that in PUPer_xlate)
  Tchar=0,Tshort, Tint, Tlong, Tlonglong,
  Tuchar,Tushort,Tuint,Tulong, Tulonglong,
  Tfloat,Tdouble,Tlongdouble,
  Tbool,
  Tbyte,
  Tsync,
  dataType_last //<- for setting table lengths, etc.
} dataType;

//This should be a 1-byte unsigned type
typedef unsigned char myByte;

//Forward declarations
class er;
class able;
class xlater;

//Used for out-of-order unpacking
class seekBlock {
	enum {maxSections=3};
	int secTab[maxSections+1];//The start of each seek section
	int nSec;//Number of sections; current section #
	int secTabOff;//Start of the section table, relative to the seek block
	er &p;
	CmiBool hasEnded;
public:
	//Constructor
	seekBlock(er &Np,int nSections);
	//Destructor
	~seekBlock();

	//Seek to the given section number (0-based, less than nSections)
	void seek(int toSection);
	//Finish with this seeker (must be called)
	void endBlock(void);

	//An evil hack to avoid inheritance and virtual functions among seekers--
	// stores the PUP::er specific block start information.
	union {
		int off;
		long loff;
		const myByte *cptr;
		myByte *ptr;
		void *vptr;
	} data;
};

//The abstract base class:  PUP::er.
class er {
 private:
  er(const er &p);//You don't want to copy PUP::er's.
 protected:
  enum {IS_DELETING =0x0008, IS_USERLEVEL =0x0004};
  enum {IS_SIZING   =0x0100,
  	IS_PACKING  =0x0200,
        IS_UNPACKING=0x0400,
        TYPE_MASK   =0xFF00};
  unsigned int PUP_er_state;
#if CMK_EXPLICIT
  explicit /* Makes constructor below behave better */
#endif
           er(unsigned int inType) //You don't want to create raw PUP::er's.
              {PUP_er_state=inType;}
 public:
  virtual ~er();//<- does nothing, but might be needed by some child

  //State queries (exactly one of these will be true)
  CmiBool isSizing(void) const {return (PUP_er_state&IS_SIZING)!=0?CmiTrue:CmiFalse;}
  CmiBool isPacking(void) const {return (PUP_er_state&IS_PACKING)!=0?CmiTrue:CmiFalse;}
  CmiBool isUnpacking(void) const {return (PUP_er_state&IS_UNPACKING)!=0?CmiTrue:CmiFalse;}
  char *  typeString() const;

  //This indicates that the pup routine should free memory during packing.
  void becomeDeleting(void) {PUP_er_state|=IS_DELETING;}
  CmiBool isDeleting(void) const {return (PUP_er_state&IS_DELETING)!=0?CmiTrue:CmiFalse;}

  //This indicates that the pup routine should not call system objects' pups.
  void becomeUserlevel(void) {PUP_er_state|=IS_USERLEVEL;}
  CmiBool isUserlevel(void) const {return (PUP_er_state&IS_USERLEVEL)!=0?CmiTrue:CmiFalse;}

//For single elements, pretend it's an array containing one element
  void operator()(signed char &v)     {(*this)(&v,1);}
#if CMK_SIGNEDCHAR_DIFF_CHAR
  void operator()(char &v)            {(*this)(&v,1);}
#endif
  void operator()(short &v)           {(*this)(&v,1);}
  void operator()(int &v)             {(*this)(&v,1);}
  void operator()(long &v)            {(*this)(&v,1);}
  void operator()(unsigned char &v)   {(*this)(&v,1);}
  void operator()(unsigned short &v)  {(*this)(&v,1);}
  void operator()(unsigned int &v)    {(*this)(&v,1);}
  void operator()(unsigned long &v)   {(*this)(&v,1);}
  void operator()(float &v)           {(*this)(&v,1);}
  void operator()(double &v)          {(*this)(&v,1);}
#if CMK_LONG_DOUBLE_DEFINED
  void operator()(long double &v)     {(*this)(&v,1);}
#endif
  void operator()(CmiBool &v)         {(*this)(&v,1);}
#ifdef CMK_PUP_LONG_LONG
  void operator()(CMK_PUP_LONG_LONG &v) {(*this)(&v,1);}
  void operator()(unsigned CMK_PUP_LONG_LONG &v) {(*this)(&v,1);}
#endif

//For arrays:
  //Integral types:
  void operator()(signed char *a,int nItems)
    {bytes((void *)a,nItems,sizeof(signed char),Tchar);}
#if CMK_SIGNEDCHAR_DIFF_CHAR
  void operator()(char *a,int nItems)
    {bytes((void *)a,nItems,sizeof(char),Tchar);}
#endif
  void operator()(short *a,int nItems)
    {bytes((void *)a,nItems,sizeof(short),Tshort);}
  void operator()(int *a,int nItems)
    {bytes((void *)a,nItems,sizeof(int),Tint);}
  void operator()(long *a,int nItems)
    {bytes((void *)a,nItems,sizeof(long),Tlong);}

  //Unsigned integral types:
  void operator()(unsigned char *a,int nItems)
    {bytes((void *)a,nItems,sizeof(unsigned char),Tuchar);}
  void operator()(unsigned short *a,int nItems)
    {bytes((void *)a,nItems,sizeof(unsigned short),Tushort);}
  void operator()(unsigned int *a,int nItems)
    {bytes((void *)a,nItems,sizeof(unsigned int),Tuint);}
  void operator()(unsigned long *a,int nItems)
    {bytes((void *)a,nItems,sizeof(unsigned long),Tulong);}

  //Floating-point types:
  void operator()(float *a,int nItems)
    {bytes((void *)a,nItems,sizeof(float),Tfloat);}
  void operator()(double *a,int nItems)
    {bytes((void *)a,nItems,sizeof(double),Tdouble);}

#if CMK_LONG_DOUBLE_DEFINED
  void operator()(long double *a,int nItems)
    {bytes((void *)a,nItems,sizeof(long double),Tlongdouble);}
#endif

  //For bools:
  void operator()(CmiBool *a,int nItems)
    {bytes((void *)a,nItems,sizeof(CmiBool),Tbool);}

#ifdef CMK_PUP_LONG_LONG
  void operator()(CMK_PUP_LONG_LONG *a,int nItems)
    {bytes((void *)a,nItems,sizeof(CMK_PUP_LONG_LONG),Tlonglong);}
  void operator()(unsigned CMK_PUP_LONG_LONG *a,int nItems)
    {bytes((void *)a,nItems,sizeof(unsigned CMK_PUP_LONG_LONG),Tulonglong);}
#endif

  //For raw memory (n gives number of bytes)
  void operator()(void *a,int nBytes)
    {bytes((void *)a,nBytes,1,Tbyte);}

  //For allocatable objects (system will new/delete object and call pup routine)
  void operator()(able** a)
    {object(a);}
  //For pre- or stack-allocated PUP::able objects-- just call object's pup
  void operator()(able& a);

  //A descriptive (but entirely optional) human-readable comment field
  virtual void comment(const char *message);

  //A 32-bit synchronization marker (not human readable)
  virtual void synchronize(unsigned int m);

 protected:
  //Generic bottleneck: pack/unpack n items of size itemSize
  // and data type t from p.  Desc describes the data item
  friend class xlater;
  virtual void bytes(void *p,int n,size_t itemSize,dataType t) =0;
  virtual void object(able** a);

  //For seeking (pack/unpack in different orders)
  friend class seekBlock;
  virtual void impl_startSeek(seekBlock &s); /*Begin a seeking block*/
  virtual int impl_tell(seekBlock &s); /*Give the current offset*/
  virtual void impl_seek(seekBlock &s,int off); /*Seek to the given offset*/
  virtual void impl_endSeek(seekBlock &s);/*End a seeking block*/
};

/************** PUP::er -- Sizer ******************/
//For finding the number of bytes to pack (e.g., to preallocate a memory buffer)
class sizer : public er {
 protected:
  int nBytes;
  //Generic bottleneck: n items of size itemSize
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Write data to the given buffer
  sizer(void):er(IS_SIZING) {nBytes=0;}
  
  //Return the current number of bytes to be packed
  int size(void) const {return nBytes;}
};

template <class T>
inline int size(T &t) {
	PUP::sizer p; p|t; return p.size();
}

/********** PUP::er -- Binary memory buffer pack/unpack *********/
class mem : public er { //Memory-buffer packers and unpackers
 protected:
  myByte *origBuf;//Start of memory buffer
  myByte *buf;//Memory buffer (stuff gets packed into/out of here)
  mem(unsigned int type,myByte *Nbuf):er(type),origBuf(Nbuf),buf(Nbuf) {}

  //For seeking (pack/unpack in different orders)
  virtual void impl_startSeek(seekBlock &s); /*Begin a seeking block*/
  virtual int impl_tell(seekBlock &s); /*Give the current offset*/
  virtual void impl_seek(seekBlock &s,int off); /*Seek to the given offset*/
 public:
  //Return the current number of buffer bytes used
  int size(void) const {return buf-origBuf;}
};

//For packing into a preallocated, presized memory buffer
class toMem : public mem {
 protected:
  //Generic bottleneck: pack n items of size itemSize from p.
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Write data to the given buffer
  toMem(void *Nbuf):mem(IS_PACKING,(myByte *)Nbuf) {}
};
template <class T>
inline void toMemBuf(T &t,void *buf,int len) {
	PUP::toMem p(buf);
	p|t;
	if (p.size()!=len) CmiAbort("Size mismatch during PUP::toMemBuf!\n"
		"This means your pup routine doesn't match during sizing and packing");
}

//For unpacking from a memory buffer
class fromMem : public mem {
 protected:
  //Generic bottleneck: unpack n items of size itemSize from p.
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Read data from the given buffer
  fromMem(const void *Nbuf):mem(IS_UNPACKING,(myByte *)Nbuf) {}
};
template <class T>
inline void fromMemBuf(T &t,void *buf,int len) {
	PUP::fromMem p(buf);
	p|t;
	if (p.size()!=len) CmiAbort("Size mismatch during PUP::fromMemBuf!\n"
		"This means your pup routine doesn't match during packing and unpacking");
}

/********** PUP::er -- Binary disk file pack/unpack *********/
class disk : public er {
 protected:
  FILE *F;//Disk file to read from/write to
  disk(unsigned int type,FILE *f):er(type),F(f) {}

  //For seeking (pack/unpack in different orders)
  virtual void impl_startSeek(seekBlock &s); /*Begin a seeking block*/
  virtual int impl_tell(seekBlock &s); /*Give the current offset*/
  virtual void impl_seek(seekBlock &s,int off); /*Seek to the given offset*/
};

//For packing to a disk file
class toDisk : public disk {
 protected:
  //Generic bottleneck: pack n items of size itemSize from p.
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Write data to the given file pointer
  // (must be opened for binary write)
  // You must close the file yourself when done.
  toDisk(FILE *f):disk(IS_PACKING,f) {}
};

//For unpacking from a disk file
class fromDisk : public disk {
 protected:
  //Generic bottleneck: unpack n items of size itemSize from p.
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Write data to the given file pointer 
  // (must be opened for binary read)
  // You must close the file yourself when done.
  fromDisk(FILE *f):disk(IS_UNPACKING,f) {}
};



/************** PUP::er -- Text *****************/
class toTextUtil : public er {
 private:
  char *cur; /*Current output buffer*/
  int level; /*Indentation distance*/
  void beginEnv(const char *type,int n=0);
  void endEnv(const char *type);
  char *beginLine(void);
  void endLine(void);
 protected:
  virtual char *advance(char *cur)=0; /*Consume current buffer and return next*/
  toTextUtil(unsigned int inType,char *buf);
 public:
  virtual void comment(const char *message);
  virtual void synchronize(unsigned int m);
 protected:
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
  virtual void object(able** a);
};
/* Return the number of characters, including terminating NULL */
class sizerText : public toTextUtil {
 private:
  char line[1000];
  int charCount; /*Total characters seen so far (not including NULL) */
 protected:
  virtual char *advance(char *cur);
 public:
  sizerText(void);
  int size(void) const {return charCount+1; /*add NULL*/ }
};
/* Copy data to this C string, including terminating NULL. */
class toText : public toTextUtil {
 private:
  char *buf;
  int charCount; /*Total characters written so far (not including NULL) */
 protected:
  virtual char *advance(char *cur);
 public:
  toText(char *outStr);
  int size(void) const {return charCount+1; /*add NULL*/ }
};

class toTextFile : public er {
 protected:
  FILE *f;
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  //Begin writing to this file, which should be opened for ascii write.
  // You must close the file yourself when done.
  toTextFile(FILE *f_) :er(IS_PACKING), f(f_) {}
  virtual void comment(const char *message);
};
class fromTextFile : public er {
 protected:
  FILE *f;
  int readInt(const char *fmt="%d");
  unsigned int readUint(const char *fmt="%u");
  double readDouble(void);
  
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
  virtual void parseError(const char *what);
 public:
  //Begin writing to this file, which should be opened for ascii read.
  // You must close the file yourself when done.
  fromTextFile(FILE *f_) :er(IS_UNPACKING), f(f_) {}
  virtual void comment(const char *message);
};

/********** PUP::er -- Heterogenous machine pack/unpack *********/
//This object describes the data representation of a machine.
class machineInfo {
 public:
  typedef unsigned char myByte;
  myByte magic[4];//Magic number (to identify machineInfo structs)
  myByte version;//0-- current version

  myByte intBytes[4]; //<- sizeof(char,short,int,long)
  myByte intFormat;//0-- big endian.  1-- little endian.

  myByte floatBytes; //<- sizeof(...)
  myByte doubleBytes;
  myByte floatFormat;//0-- big endian IEEE.  1-- little endian IEEE.

  myByte boolBytes;

  myByte padding[2];//Padding to 16 bytes

  //Return true if our magic number is valid.
  CmiBool valid(void) const;
  //Return true if we differ from the current (running) machine.
  CmiBool needsConversion(void) const;
  
  //Get a machineInfo for the current machine
  static const machineInfo &current(void);
};

//For translating some odd disk/memory representation into the 
// current machine representation.  (We really only need to
// translate during unpack-- "reader makes right".)
class xlater : public er {
 protected:
  typedef void (*dataConverterFn)(int N,const myByte *in,myByte *out,int nElem);
  
  //This table is indexed by dataType, and contains an appropriate
  // conversion function to unpack a n-item array of the corresponding 
  // data type (possibly in-place).
  dataConverterFn convertFn[dataType_last];
  //Maps dataType to source machine's dataSize
  size_t convertSize[dataType_last];
  void setConverterInt(const machineInfo &m,const machineInfo &cur,
    int isUnsigned,int intType,dataType dest);
  
  er &myUnpacker;//Raw data unpacker
  //Generic bottleneck: unpack n items of size itemSize from p.
  virtual void bytes(void *p,int n,size_t itemSize,dataType t);
 public:
  xlater(const machineInfo &fromMachine, er &fromData);
 protected:
  //For seeking (pack/unpack in different orders)
  friend class seekBlock;
  virtual void impl_startSeek(seekBlock &s); /*Begin a seeking block*/
  virtual int impl_tell(seekBlock &s); /*Give the current offset*/
  virtual void impl_seek(seekBlock &s,int off); /*Seek to the given offset*/
  virtual void impl_endSeek(seekBlock &s);/*End a seeking block*/
};

/*************** PUP::able support ***************/
//The base class of system-allocatable objects with pup routines
class able {
public:
	//A globally-unique, persistent identifier for an allocatable object
	class PUP_ID {
	public:
		enum {len=8};
		unsigned char hash[len];
		PUP_ID() {}
		PUP_ID(int val) {for (int i=0;i<len;i++) hash[i]=val;}
		PUP_ID(const char *name) {setName(name);}
		void setName(const char *name);//Write name into hash
		CmiBool operator==(const PUP_ID &other) const {
			for (int i=0;i<len;i++)
				if (hash[i]!=other.hash[i])
					return CmiFalse;
			return CmiTrue;
		}
		void pup(er &p) {
			 p((void *)hash,sizeof(unsigned char)*len);
		}
		void pup(er &p) const {
			 p((void *)hash,sizeof(unsigned char)*len);
		}
	};

protected:
	able() {}
	able(CkMigrateMessage *) {}
	virtual ~able();//Virtual destructor may be needed by some child

public:
//Constructor function registration:
	typedef able* (*constructor_function)(void);
	static PUP_ID register_constructor(const char *className,
		constructor_function fn);
	static constructor_function get_constructor(const PUP_ID &id);
	virtual /*PUP::*/able *clone(void) const;

//Target methods:
	virtual void pup(er &p);
	virtual const PUP_ID &get_PUP_ID(void) const=0;

};

//Declarations to include in a PUP::able's body
#define PUPable_decl(className) \
private: \
    static PUP::able *call_PUP_constructor(void); \
    static PUP::able::PUP_ID my_PUP_ID;\
public:\
    virtual const PUP::able::PUP_ID &get_PUP_ID(void) const; \
    static void register_PUP_ID(void); \
    friend inline void operator|(PUP::er &p,className &a) {a.pup(p);}\
    friend inline void operator|(PUP::er &p,className* &a) {\
	PUP::able *pa=a;  p(&pa);  a=(className *)pa;\
    }

//Declarations to include in an abstract PUP::able's body.
//  Abstract PUP::ables do not need def or reg.
#define PUPable_abstract(className) \
public:\
    virtual const PUP::able::PUP_ID &get_PUP_ID(void) const =0; \
    friend inline void operator|(PUP::er &p,className &a) {a.pup(p);}\
    friend inline void operator|(PUP::er &p,className* &a) {\
	PUP::able *pa=a;  p(&pa);  a=(className *)pa;\
    }

//Definitions to include exactly once at file scope
#define PUPable_def(className) \
	PUP::able *className::call_PUP_constructor(void) \
		{ return new className((CkMigrateMessage *)0);}\
	const PUP::able::PUP_ID &className::get_PUP_ID(void) const\
		{ return className::my_PUP_ID; }\
	PUP::able::PUP_ID className::my_PUP_ID;\
	void className::register_PUP_ID(void)\
		{my_PUP_ID=register_constructor(#className,\
		              className::call_PUP_constructor);}\

//Code to execute exactly once at program start time
#define PUPable_reg(className) \
    className::register_PUP_ID();


};//<- End namespace PUP

inline void operator|(PUP::er &p,PUP::able &a) {a.pup(p);}
inline void operator|(PUP::er &p,PUP::able* &a) {p(&a);}

//Holds a pointer to a (possibly dynamically allocated) PUP::able.
//  Extracting the pointer hands the deletion responsibility over.
//  This is used by parameter marshalling, which doesn't work well 
//  with bare pointers.
//   CkPointer<T> t   is the parameter-marshalling equivalent of   T *t
extern "C" void CmiAbort(const char *msg);
template <class T>
class CkPointer {
	T *allocated; //Pointer that PUP dynamically allocated for us (recv only)
	T *ptr; //Read-only pointer

#if 0 /* Private (do-not-use) copy constructor.  This prevents allocated from being
         deleted twice--once in the original, and again in the copy.*/
	CkPointer(const CkPointer<T> &src); // Don't use this!
#else /* Some compilers, like gcc3, have a hideous bug that causes them to *demand*
         a public copy constructor when a class is used to initialize a const-reference
	 from a temporary.  The public copy constructor should never be called, though. */
public:
	CkPointer(const CkPointer<T> &src) {
		CmiAbort("PUPable_marshall's cannot be passed by value.  Pass them only by reference!");
	}
#endif
protected:
	T *peek(void) {return ptr;}
public:
	/// Used on the send side, and does *not* delete the object.
	CkPointer(T *src)  ///< Marshall this object.
	{ 
		allocated=0; //Don't ever delete src
		ptr=src;
	}
	
	/// Begin completely empty: used on marshalling recv side.
	CkPointer(void) { 
		ptr=allocated=0;
	}
	
	~CkPointer() { if (allocated) delete allocated; }
	
	/// Extract the object held by this class.  
	///  Deleting the pointer is now the user's responsibility
	inline operator T* () { allocated=0; return ptr; }
	
	inline void pup(PUP::er &p) {
		bool ptrWasNull=(ptr==0);
		
		PUP::able *ptr_able=ptr; // T must inherit from PUP::able!
		p|ptr_able; //Pack as a PUP::able *
		ptr=(T *)ptr_able;
		
		if (ptrWasNull) 
		{ //PUP just allocated a new object for us-- 
		  // make sure it gets deleted eventually.
			allocated=ptr;
		}
	}
	friend inline void operator|(PUP::er &p,CkPointer<T> &v) {v.pup(p);}
};
#define PUPable_marshall CkPointer

//Like CkPointer, but keeps deletion responsibility forever.
//   CkReference<T> t  is the parameter-marshalling equivalent of   T &t
template<class T>
class CkReference : private CkPointer<T> {
public:
	/// Used on the send side, and does *not* delete the object.
	CkReference(T &src)   ///< Marshall this object.
		:CkPointer<T>(&src) { }
	
	/// Begin completely empty: used on the recv side.
	CkReference(void) {}
	
	/// Look at the object held by this class.  Does *not* hand over
	/// deletion responsiblity.
	inline operator T& () { return *peek(); }
	
	inline void pup(PUP::er &p) {CkPointer<T>::pup(p);}
	
	friend inline void operator|(PUP::er &p,CkReference<T> &v) {v.pup(p);}
};

// For people that forget ::'s:
typedef PUP::er PUPer;
typedef PUP::able PUPable;

/******** PUP via pipe: another way to access PUP::ers *****
The parameter marshalling system pups each variable v using just:
     p|v;
Thus we need a "void operator|(PUP::er &p,T &v)" for all types
that work with parameter marshalling.  This operator| is often
defined by the PUPmarshall macro, below.
*/

//These versions map p|t to p(t) for all built-in types
inline void operator|(PUP::er &p,signed char &t) {p(t);}
#if CMK_SIGNEDCHAR_DIFF_CHAR
inline void operator|(PUP::er &p,char &t) {p(t);}
#endif
inline void operator|(PUP::er &p,unsigned char &t) {p(t);}
inline void operator|(PUP::er &p,short &t) {p(t);}
inline void operator|(PUP::er &p,int &t) {p(t);}
inline void operator|(PUP::er &p,long &t) {p(t);}
inline void operator|(PUP::er &p,unsigned short &t) {p(t);}
inline void operator|(PUP::er &p,unsigned int &t) {p(t);}
inline void operator|(PUP::er &p,unsigned long &t) {p(t);}
inline void operator|(PUP::er &p,float &t) {p(t);}
inline void operator|(PUP::er &p,double &t) {p(t);}
inline void operator|(PUP::er &p,CmiBool &t) {p(t);}
#if CMK_LONG_DOUBLE_DEFINED
inline void operator|(PUP::er &p,long double &t) {p(t);}
#endif
#ifdef CMK_PUP_LONG_LONG
inline void operator|(PUP::er &p,CMK_PUP_LONG_LONG &t) {p(t);}
inline void operator|(PUP::er &p,unsigned CMK_PUP_LONG_LONG &t) {p(t);}
#endif

#ifdef CK_DEFAULT_BITWISE_PUP
/// This defines "p|t" as a byte-by-byte copy, for any
///  user-defined type T that does not have a normal operator| defined.
/// It is the old, "convenient but error-prone" definition.
template <class T>
inline void operator|(PUP::er &p,T &t)
{
         p((void *)&t,sizeof(T));
}

//Map operator| to this classes' ordinary pup routine.
#define PUPmarshall(type) \
  inline void operator|(PUP::er &p,type &t) {t.pup(p);}

#else /* The usual case: !CK_DEFAULT_BITWISE_PUP */

/// This defines "p|t" to call t's pup routine, if no
///  existing operator| is found.
template <class T>
inline void operator|(PUP::er &p,T &t) {t.pup(p);}

#define PUPmarshall(type) /* empty, for backward compatability */

#endif

#define PUPmarshal(type) PUPmarshall(type) /*Support this common misspelling*/


/// Copy these classes as raw memory
#define PUPbytes(type) \
  inline void operator|(PUP::er &p,type &t) { p((void *)&t,sizeof(type)); }

#define PUPmarshallBytes(type) PUPbytes(type)


#endif //def __CK_PUP_H


