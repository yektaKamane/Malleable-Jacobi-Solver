CMK_CFLAGS="-xarch=v9a -KPIC"
CMK_CXXFLAGS="-xarch=v9a -KPIC"
CMK_CPP_CHARM='/usr/ccs/lib/cpp '
CMK_CPP_C='cc -E '
CMK_CC="cc $CMK_CFLAGS"
CMK_CC_RELIABLE="cc $CMK_CFLAGS"
CMK_CC_FASTEST="cc $CMK_CFLAGS"
CMK_CXX="CC -library=Cstd $CMK_CXXFLAGS"
CMK_CXXPP='CC -E '
CMK_CF77='f77'
CMK_CF90='f90 -C -stackvar -g -xarch=v9a -KPIC'
CMK_C_DEBUG='-g'
CMK_C_OPTIMIZE='-xO5'
CMK_CXX_DEBUG='-g'
CMK_CXX_OPTIMIZE='-xO5'
CMK_LD="cc $CMK_CFLAGS "
CMK_LDXX="CC -library=Cstd $CMK_CXXFLAGS"
CMK_AR='CC -xar -o'
CMK_LIBS='-lnsl -lsocket -lckqt'
CMK_LD_SHARED="-G"
CMK_SEQ_LIBS='-lnsl -lsocket'
CMK_SEQ_CC="cc $CMK_CFLAGS "
CMK_SEQ_LD="cc $CMK_CFLAGS "
CMK_SEQ_CXX="CC -library=Cstd $CMK_CXXFLAGS "
CMK_SEQ_LDXX="CC -library=Cstd $CMK_CXXFLAGS "
CMK_CPP_SUFFIX="cc"
CMK_XLATPP='charmxlat++ -w -p '
CMK_QT='solaris-cc64'
CMK_F90LIBS='-lfsu -lsunmath -lfsumai -lfminvai -lfmaxvai -lfui -lfai'
CMK_MOD_EXT="mod"
CMK_F90_USE_MODDIR=1
