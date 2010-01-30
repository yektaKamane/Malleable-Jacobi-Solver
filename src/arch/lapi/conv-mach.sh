COMMENT="Enable 32-bit mode (-q32)"
if [ "$OBJECT_MODE" != "32" ]
then
        echo "WARNING: Charm built in 32-bit mode, but OBJECT_MODE=$OBJECT_MODE"
fi

CMK_CPP_CHARM='/usr/lib/cpp -P -D_NO_PROTO '
CMK_CPP_C='/usr/lib/cpp -P -D_NO_PROTO '
CMK_LDRO='ld -r -o '
CMK_LDRO_WORKS=0
CMK_CC='mpcc_r -q32 -qcpluscmt '
CMK_CXX='mpCC_r -q32 -qstaticinline '
CMK_CXXPP='xlC -q32 -E '
CMK_LD="mpcc_r -q32 -brtl "
CMK_LDXX="mpCC_r -q32 -brtl "
CMK_CF77='mpxlf_r -q32 '
CMK_CF90='mpxlf90_r -q32 -qsuffix=f=f90'
CMK_CF90_FIXED='mpxlf90_r -q32 '
CMK_C_OPTIMIZE='-O3 -qstrict -Q  '
CMK_CXX_OPTIMIZE='-O3 -qstrict -Q '
CMK_AR='ar -X 32 cq'
CMK_RANLIB='true'
CMK_LIBS="-lckqt -lhC -llapi_r"
CMK_LD_SHARED='-G'
CMK_SEQ_LIBS=''
CMK_SEQ_CC='xlc_r -q32 '
CMK_SEQ_LD='xlc_r -q32 '
CMK_SEQ_CXX='xlC_r -q32 -qstaticinline '
CMK_SEQ_LDXX='xlC_r -q32 '
CMK_NATIVE_CC='xlc_r -q32 '
CMK_NATIVE_LD='xlc_r -q32 '
CMK_NATIVE_CXX='xlC_r -q32 -qstaticinline '
CMK_NATIVE_LDXX='xlC_r -q32 '
CMK_NM='/bin/nm -X 32'
CMK_NM_FILTER="grep ^_CK_ | cut -f 1 -d ' '"
CMK_QT='aix'
CMK_XIOPTS=''
CMK_MOD_EXT='mod'
CMK_F90LIBS='-lxlf90_r -lhC -llapi_r'
