#cmake_minimum_required(VERSION 2.8)
#project(METIS)


# Add flags known to be needed for library. Assume that DEBUG has potentially lower optimization
add_compile_options(-O2)

# Create libparmetis
add_library(metis  coarsen.c fm.c initpart.c match.c ccgraph.c
       pmetis.c pqueue.c refine.c util.c timing.c debug.c
       bucketsort.c graph.c stat.c kmetis.c kwayrefine.c
       kwayfm.c balance.c ometis.c srefine.c sfm.c separator.c
       mincover.c mmd.c mesh.c meshpart.c frename.c fortran.c
       myqsort.c compress.c parmetis.c estmem.c
       mpmetis.c mcoarsen.c mmatch.c minitpart.c mbalance.c
       mutil.c mkmetis.c mkwayrefine.c mkwayfmh.c
       mrefine2.c minitpart2.c mbalance2.c mfm2.c
       kvmetis.c kwayvolrefine.c kwayvolfm.c subdomains.c
       mfm.c memory.c mrefine.c checkgraph.c)



target_compile_definitions(metis PRIVATE -D_Float128=__float128)


# Add include directories.
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})


