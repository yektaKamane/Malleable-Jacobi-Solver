CMK_USE_PAPI=true
USE_SPP_PAPI=true
#you should run module load papi
if test "x$PAT_BUILD_PAPI_BASEDIR" != "x"
then
PAPI_LIBDIR="$PAT_BUILD_PAPI_BASEDIR/lib"
PAPI_INCDIR="$PAT_BUILD_PAPI_BASEDIR/include"
else
PAPI_LIBDIR="/opt/cray/papi/4.3.0.1/perf_events/no-cuda/lib"
PAPI_INCDIR="/opt/cray/papi/4.3.0.1/perf_events/no-cuda/include"
fi
CMK_INCDIR="$CMK_INCDIR -I$PAPI_INCDIR"
CMK_LIBDIR="-L $PAPI_LIBDIR"
CMK_LD="$CMK_LD -Wl,-rpath,$PAPI_LIBDIR"
CMK_LDXX="$CMK_LDXX -Wl,-rpath,$PAPI_LIBDIR" 
CMK_LIBS="$CMK_LIBS -lpapi"
