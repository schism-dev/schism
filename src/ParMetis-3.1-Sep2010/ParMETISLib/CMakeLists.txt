# Include directories for library code.
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/..)

# Add flags known to be needed for library. Assume that DEBUG has potentially lower optimization
add_compile_options(-O2)

# Create libparmetis
add_library(parmetis comm.c util.c debug.c setup.c grsetup.c timer.c
       node_refine.c initmsection.c order.c
       xyzpart.c pspases.c frename.c
       iintsort.c iidxsort.c ikeysort.c ikeyvalsort.c
       kmetis.c gkmetis.c ometis.c
       initpart.c match.c coarsen.c
       kwayfm.c kwayrefine.c kwaybalance.c
       remap.c stat.c fpqueue.c
       ametis.c rmetis.c lmatch.c initbalance.c
       mdiffusion.c diffutil.c wave.c
       csrmatch.c redomylink.c balancemylink.c
       selectq.c akwayfm.c serial.c move.c
       mmetis.c mesh.c memory.c weird.c backcompat.c)



#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLATE_SUB_CXX_FLAG")


# Link with metis and MPI libraries.
target_link_libraries(parmetis metis ${MPI_C_LIBRARIES})
set_target_properties(parmetis PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
target_compile_definitions(parmetis PRIVATE -D_Float128=__float128)
