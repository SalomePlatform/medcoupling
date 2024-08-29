find_package(cppunit REQUIRED)

add_subdirectory(interp_kernel)
add_subdirectory(medcoupling_core)

if(NOT MEDCOUPLING_MICROMED)
  add_subdirectory(medfile_api)
endif()

if(MEDCOUPLING_USE_MPI)
  add_subdirectory(para_medmem)
endif(MEDCOUPLING_USE_MPI)
