CMK_CPP_CHARM="/lib/cpp -P"
CMK_CPP_C="gcc -E"
CMK_CC="gcc "
CMK_CXX="g++ "
CMK_CXXPP="$CMK_CXX -x c++ -E "
CMK_RANLIB="ranlib"
CMK_LIBS="-lckqt"
#CMK_LD="$CMK_CC -Wl,--allow-multiple-definition "
#CMK_LDXX="$CMK_CXX -Wl,--allow-multiple-definition "
CMK_LD="$CMK_CC "
CMK_LDXX="$CMK_CXX "
CMK_LD_SHARED="-shared"
CMK_LD_LIBRARY_PATH="-Wl,-rpath,$CHARMLIBSO/"
CMK_XIOPTS=""

# fortran compiler Absoft or gnu f95
CMK_CF77="g77 "
CMK_F77LIBS="-lg2c "
CMK_CF90=`which f90 2>/dev/null`
if test -n "$CMK_CF90"
then
# absoft
CMK_CF90_FIXED="$CMK_CF90 -W132 "
CMK_F90LIBS="-L/usr/absoft/lib -L/opt/absoft/lib -lf90math -lfio -lU77 -lf77math "
CMK_MOD_NAME_ALLCAPS=1
CMK_MOD_EXT="mod"
CMK_F90_USE_MODDIR=1
CMK_F90_MODINC="-p"
else
# gnu f95
  CMK_CF90=`which f95 2>/dev/null`
  if test -n "$CMK_CF90"
  then
    CMK_FPP="/lib/cpp -P -CC"
    CMK_CF90_FIXED="$CMK_CF90 -ffixed-form "
    CMK_F90LIBS="-lgfortran "
    CMK_MOD_NAME_ALLCAPS=
    CMK_MOD_EXT="mod"
    CMK_F90_USE_MODDIR=1
    CMK_F90_MODINC="-I"
  fi
fi
