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

macro(tohexa num result)
  set(_hexa_map a b c d e f)
  if(${num} LESS 10)
    set(${result} "0${num}")
  else()
    math(EXPR _res "${num}-10" )
    list(GET _hexa_map ${_res} _out)
    set(${result} "0${_out}")
  endif()
  message("Transforming " ${num} " into hex " ${${result}})
endmacro()

tohexa(${${PROJECT_NAME}_VERSION_MAJOR} _major)
tohexa(${${PROJECT_NAME}_VERSION_MINOR} _minor)
tohexa(${${PROJECT_NAME}_VERSION_PATCH} _patch)
set(${PROJECT_NAME}_XVERSION "0x${_major}${_minor}${_patch}")

configure_file(include/MEDCoupling_version.h.in include/MEDCoupling_version.h)

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
