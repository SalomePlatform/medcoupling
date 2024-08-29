if(MEDCOUPLING_ENABLE_RENUMBER)
  add_subdirectory(renumberer)
endif()

if(MEDCOUPLING_ENABLE_PARTITIONER)
  add_subdirectory(partitioner)
endif()

if(MEDCOUPLING_ENABLE_STAT)
  add_subdirectory(stat)
endif()
