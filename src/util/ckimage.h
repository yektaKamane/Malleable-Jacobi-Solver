/*
2D flat image class:
This class represents a 2D raster image; a rectangular 2D
array of pixels.

Orion Sky Lawlor, olawlor@acm.org, 5/15/2002
*/
#ifndef __CK_IMAGE_H
#define __CK_IMAGE_H

#include "pup.h"

#undef min
#undef max
inline int min(int a,int b) {return (a<b)?a:b;}
inline int max(int a,int b) {return (a>b)?a:b;}
class CkRect {
public:
	int l,r; //X boundaries of rectangle
	int t,b; //Y boundaries of rectangle
	CkRect() {l=r=t=b=-1;}
	CkRect(int l_,int t_,int r_,int b_) 
		:l(l_), r(r_), t(t_), b(b_) {}
	CkRect(int w,int h)
		:l(0), r(w), t(0), b(h) {}
	//Default copy constructor, assignment operator
	int wid(void) const {return r-l;}
	int ht(void) const {return b-t;}
	int getWidth(void) const {return r-l;}
	int getHeight(void) const {return b-t;}
	inline int operator==(const CkRect &a) 
		{return l==a.l && r==a.r && t==a.t && b==a.b;}
	CkRect getUnion(const CkRect &a) {
		return CkRect(min(l,a.l),min(t,a.t), max(r,a.r),max(b,a.b));
	}
	CkRect getIntersect(const CkRect &a) {
		return CkRect(max(l,a.l),max(t,a.t), min(r,a.r),min(b,a.b));
	}
	CkRect getShift(int dx,int dy) {
		return CkRect(l+dx,t+dy,r+dx,b+dy);
	}
	CmiBool isEmpty(void) const {return (CmiBool)((l>=r) || (t>=b));}
	CmiBool inbounds(int x,int y) const {
		if (x<l || x>=r) return CmiFalse;
		if (y<t || y>=b) return CmiFalse;
		return CmiTrue;
	}
	void makeEmpty(void) {l=t=1000000000; b=r=-1000000000;}
	void empty(void) {makeEmpty();}
	void add(int x,int y) {
		l=min(x,l); r=max(x,r);
		t=min(y,t); b=max(y,b);
	}
	void enlarge(int dx,int dy) {
		l-=dx; r+=dx; t-=dy; b+=dy;
	}
	void zero(void) {l=r=t=b=0;}
	int area(void) const {return (r-l)*(b-t);}
	
	void pup(PUP::er &p) {
		p|l; p|r; p|t; p|b;
	}
};
PUPmarshall(CkRect);

/*This class describes a flat byte array interpreted as an image.
Pixels are stored first by color (r,g,b), then by row.  
*/
class CkImage {
public:
	//This is actually the data type of a color channel, not a pixel
	typedef unsigned char pixel_t;
private:
	int row,colors; //pixel_ts per line, pixel_ts per pixel
	int wid,ht; //Image size: cols and rows
	pixel_t *data; //Image pixel data
	
	CkImage(const CkImage &im) ; //DO NOT USE
	void operator=(const CkImage &im);
public:
	CkImage() {row=colors=wid=ht=-1; data=NULL;}
	CkImage(int w_,int h_,int colors_,pixel_t *data_)
		:row(w_*colors_), colors(colors_), wid(w_), ht(h_), data(data_) {}
	
	pixel_t *getData(void) {return data;}
	void setData(pixel_t *d) {data=d;}
	
	CkRect getRect(void) const {return CkRect(0,0,wid,ht);}
	int getColors(void) const {return colors;}
	
	int getWidth(void) const {return wid;}
	int getHeight(void) const {return ht;}
	
	//Copy the pixel at src onto the one at dest
	inline void copyPixel(const pixel_t *src,pixel_t *dest) {
		for (int i=0;i<colors;i++)
			dest[i]=src[i];
	}
	//Add the pixel at src to the one at dest, ignoring overflow
	inline void addPixel(const pixel_t *src,pixel_t *dest) {
		for (int i=0;i<colors;i++)
			dest[i]+=src[i];
	}
	//Add the pixel at src to the one at dest, clipping instead of overflowing
	inline void addPixelClip(const pixel_t *src,pixel_t *dest,
		const pixel_t *clip) 
	{
		for (int i=0;i<colors;i++)
			dest[i]=clip[(int)dest[i]+(int)src[i]];
	}
	
	
	//Get a pixel
	inline pixel_t *getPixel(int x,int y) {return data+x*colors+y*row;}
	inline const pixel_t *getPixel(int x,int y) const {return data+x*colors+y*row;}
	
	
	/*
	 Clip out this subregion of this image-- make us a subregion
	 */
	void window(const CkRect &src) {
		data+=src.t*row+src.l*colors;
		wid=src.wid(); ht=src.ht();
	}
	
	/*
	Zero out this image-- make it all black.
	*/
	void clear(void);
	
	/*
	 Copy all of src onto this image starting at (x,y).
	 */
	void put(int sx,int sy,const CkImage &src); 
	
	/*
	 Add all of src onto this image starting at (x,y).
	 */
	void add(int sx,int sy,const CkImage &src);
	/*
	 Add all of src onto this image starting at (x,y), clipping
         values instead of overflowing.
	 */
	void addClip(int sx,int sy,const CkImage &src,const pixel_t *clip);
	
	//Allocate clipping array for above routine
	static pixel_t *newClip(void);
	
	//Pup only the image *size*, not the image *data*.
	void pup(PUP::er &p) {
		p|wid; p|ht; p|colors; p|row;
	}
};
PUPmarshall(CkImage);


//A heap-allocated image
class CkAllocImage : public CkImage {
	pixel_t *allocData;
public:
	CkAllocImage() {allocData=NULL;}
	CkAllocImage(int w,int h,int c)
		:CkImage(w,h,c,new pixel_t[w*h*c]) 
	{
		allocData=getData();
	}
	~CkAllocImage() {delete[] allocData;}
	
	//Pup both image size as well as image data.
	void pup(PUP::er &p);
};
PUPmarshall(CkAllocImage);


#endif

