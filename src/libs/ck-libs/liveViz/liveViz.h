/*
  Data types, function prototypes,  etc. exported by liveViz.
  This layer does image assembly, and is the most commonly-used 
  interface to liveViz.
 */
#ifndef __UIUC_CHARM_LIVEVIZ_H
#define __UIUC_CHARM_LIVEVIZ_H

#include "liveViz0.h"
#include "ckimage.h"
#include "colorScale.h"
#include "ImageData.h"

/********************** LiveVizPoll **********************
These declarations should probably live in a header named "liveVizPoll.h"
*/
#include "liveVizPoll.decl.h"


// liveVizPollMode controls how the image is reassembled:
typedef enum {
/* This mode responds to each request immediately, ignoring
any differences between timesteps ("time skew").  This is the
fastest, simplest mode, and hence the default.
*/
        liveVizPollMode_skew=1234,

/* Explicitly synchronize timesteps, which may require several
potential images to be generated and buffered while synchronizing.
This mode is not yet implemented.
*/
        liveVizPollMode_synch
/* other modes may be added in the future */
} liveVizPollMode;

/*
Initialize the poll mode of liveViz.  This routine should
be called from main::main or just before creating the array
you wish to make liveVizPoll calls from.  For example:

        CkArrayOptions opts(nElements);
        liveVizPollInit(liveVizConfig(false,true),opts);
        CProxy_myArr arr=CProxy_myArr::ckNew(opts);

[note to implementor]
liveVizPollInit exists to allow your library to attach listener objects
to the new array.  The listeners would recieve ordinary liveViz
request callbacks, queue up the requests until the next
liveVizPoll, and perform any buffering or synchronization needed.
*/
void liveVizPollInit(const liveVizConfig &cfg,
                     CkArrayOptions &opts,
                     liveVizPollMode mode=liveVizPollMode_skew);

void liveVizPollSync(ArrayElement *from);


class liveVizPollRequestMsg : public CMessage_liveVizPollRequestMsg {
public:
  liveVizRequest3d req;
  liveVizPollRequestMsg() {}
  liveVizPollRequestMsg(const liveVizRequest3d &req_) :req(req_) {}
};


/*
Asks liveViz if there are any requests for images of this
timestep.  If there are no requests, this routine returns NULL.
There may be multiple requests for one timestep, so this routine
should be called again and again until it returns NULL.

Timesteps need not be integers, but must strictly increase.
liveVizPoll is collective, in the sense that for any timestep t,
either nobody calls liveVizPoll(this,t) or else everybody calls
liveVizPoll(this,t).


[note to implementor]
This routine can be as simple as:
        return listeners.ckLocal(from->thisIndexMax)->popRequest();
or it could initiate some complex synchronization algorithm.

Timestep is double-precision to avoid integer overflow if the
algorithm performs more than 2^32 steps (if each step takes
1us, 2^32 steps would take 1 hour).
*/
liveVizPollRequestMsg *liveVizPoll(ArrayElement *from,double timestep);


/*
Responds to a previously poll'd request.  The latter parameters
have exactly the same meaning as with liveVizDeposit.

Each non-NULL response from liveVizPoll should be followed by
a call to this routine.  An example of the proper way to
use these two routines together is:

        liveVizPollRequestMsg *req;
        while (NULL!=(req=liveVizPoll(this,timestep))) {
                //...prepare image for request...
                liveVizPollDeposit(this,timestep,req,
                        startx,starty,sizex,sizey,img);
        }

[note to implementor]
This routine may immediately call liveVizDeposit, or could buffer
the deposit until some synchronization is reached.  The exact
implementation depends strongly on which liveVizPollMode is
selected.
*/
void liveVizPollDeposit(ArrayElement *from,
                        double timestep,
			const liveVizRequest &req,
			int startx, int starty,
			int sizex, int sizey, const byte * imageData);

inline void liveVizPollDeposit(ArrayElement *from,
                        double timestep,
			liveVizPollRequestMsg *reqMsg,
			int startx, int starty,
			int sizex, int sizey, const byte * imageData)
{
  if ( ! reqMsg ) CkError("LiveViz error: User passed a null message! Bad user!");
  liveVizPollDeposit(from, timestep, reqMsg->req, startx, starty, sizex, sizey, imageData);
  delete reqMsg;
}


/********************** LiveViz ***********************/
#include "liveViz.decl.h"

/*
  Start liveViz.  This routine should be called once on processor 0
  when you are ready to begin receiving image requests.

  The arrayID is the array that will make deposits into liveViz.

  The callback is signalled each time a client requests an image.
  The image parameters are passed in as a "liveVizRequestMsg *" message.
*/
void liveVizInit(const liveVizConfig &cfg, CkArrayID a, CkCallback c);

class liveVizRequestMsg : public CMessage_liveVizRequestMsg {
public:
	liveVizRequest3d req;
	liveVizRequestMsg() {}
	liveVizRequestMsg(const liveVizRequest3d &req_) :req(req_) {}
};

/*
  Deposit a (sizex x sizey) pixel portion of the final image,
  starting at pixel (startx,starty) in the final image.
  The "client" pointer is used to perform reductions, it's
  normally "this".  Each array element must call deposit, even
  if it's just an empty deposit, like:
  	liveVizDeposit(0,0, 0,0, NULL, this);
*/
void liveVizDeposit(const liveVizRequest &req,
		    int startx, int starty,
		    int sizex, int sizey, const byte * imageData,
		    ArrayElement* client);


//As above, but taking a message instead of a request:
inline void liveVizDeposit(liveVizRequestMsg *reqMsg,
		    int startx, int starty,
		    int sizex, int sizey, const byte * imageData,
		    ArrayElement* client)
{
	liveVizDeposit(reqMsg->req,startx,starty,sizex,sizey,imageData,client);
	delete reqMsg;
}

#endif /* def(thisHeader) */
