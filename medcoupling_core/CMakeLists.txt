cmake_minimum_required(VERSION 3.24)

project(
  medcoupling_core
  VERSION 10.0.0
  DESCRIPTION
    "medcoupling_core is a C++ library defining a data format for arrays, meshes, fields and some other usefull related things. It has no dependencies and will keep being that way."
  HOMEPAGE_URL "https://github.com/SalomePlatform/medcoupling"
  LANGUAGES CXX)

option(MEDCOUPLING_ENABLE_PYTHON "Build PYTHON bindings." ON)
option(MEDCOUPLING_USE_64BIT_IDS "Size of IDs to refer cells and nodes." ON)
option(MEDCOUPLING_USE_MPI "Adds MPI capabilities to medcoupling_core" OFF)
option(MEDCOUPLING_MICROMED "Switch off med dependency" OFF)

add_subdirectory(source)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

if(MEDCOUPLING_ENABLE_PYTHON)
  find_package(
    Python
    COMPONENTS Interpreter Development.Module NumPy
    REQUIRED)

  find_package(
    SWIG
    COMPONENTS python
    REQUIRED)

  include(UseSWIG)

  list(APPEND CMAKE_SWIG_FLAGS "-DWITH_NUMPY")
  list(APPEND CMAKE_SWIG_FLAGS "-DWITH_SCIPY")

  if(MEDCOUPLING_USE_64BIT_IDS)
    list(APPEND CMAKE_SWIG_FLAGS "-DMEDCOUPLING_USE_64BIT_IDS")
  endif()

  add_subdirectory(python)
endif()
