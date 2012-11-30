###
#
# @file CMakeLists.txt
#
# CMake Information
# PLASMA is a software package provided by Univ. of Tennessee,
# Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.1.0
# @author Asim YarKhan
# @date 2010-09-20
#
###

set( quark_SRCS icl_hash.c icl_list.c quarkos.c quark.c )
set( quark_HEADERS quark.h quark_unpack_args.h icl_hash.h icl_list.h quarkwinthread.h bsd_queue.h bsd_tree.h quark_trace.h )

## This needs to be enabled if QUARK is being used outside of PLASMA
## This simple duplicates the Windows pthreads implementation in plasmawinthread.c.
#if(WIN32)
#  set( quark_SRCS ${quark_SRCS} quarkwinthread.c )
#endif( WIN32 )

add_library( quark ${quark_SRCS} )

INSTALL( TARGETS quark ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin )
INSTALL( FILES ${quark_SRCS} ${quark_HEADERS} quarkwinthread.c quarkos-hwloc.c CMakeLists.txt LICENSE DESTINATION quark )
INSTALL( FILES ${quark_HEADERS} DESTINATION quark )

