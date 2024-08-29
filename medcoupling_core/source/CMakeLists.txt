add_subdirectory(interp_kernel)
add_subdirectory(medcoupling_core)
add_subdirectory(icoco)

# medfile_api needs to be splitten in medreader/medwriter (will go into
# medcoupling_io) and polymesh => which is not dependant on Medfile
if(NOT MEDCOUPLING_MICROMED)
  add_subdirectory(medfile_api)
  if(MEDCOUPLING_USE_MPI)
    add_subdirectory(para_medfile)
  endif()
endif()

# para_medmem needs to be splitten between the para_remapper (to be moved into
# the tools) and the parallel enabled medcoupling_core classes (to be moved in
# the corresponding modules).
if(MEDCOUPLING_USE_MPI)
  add_subdirectory(para_medmem)
endif()
