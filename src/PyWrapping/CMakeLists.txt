# Copyright (C) 2017-2025  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#
# Author : Anthony Geay (EDF R&D)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

add_definitions(
  ${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS}
  -DMEDCOUPLING_GIT_SHA1=\"${MEDCOUPLING_GIT_SHA1}\"
  -DMEDCOUPLING_CONFIG_DT=\"${MEDCOUPLING_CONFIG_TIME}\")

set_source_files_properties(medcoupling.i PROPERTIES CPLUSPLUS ON)
if(WIN32)
  set_property(SOURCE medcoupling.i PROPERTY COMPILE_DEFINITIONS WIN32)
endif()
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(medcoupling.i PROPERTIES SWIG_FLAGS "-py3")
else()
  set_source_files_properties(medcoupling.i PROPERTIES SWIG_DEFINITIONS
                                                       "-shadow")
endif()
set(SWIG_MODULE_medcoupling_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
if(MEDCOUPLING_USE_64BIT_IDS)
  string(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS
         ";-DMEDCOUPLING_USE_64BIT_IDS")
endif(MEDCOUPLING_USE_64BIT_IDS)

set(medcoupling_SWIG_DPYS_FILES medcoupling.i)

include_directories(
  ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../RENUMBER
  ${CMAKE_CURRENT_SOURCE_DIR}/../RENUMBER_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDPartitioner
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDPartitioner_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../ShapeRecogn
  ${CMAKE_CURRENT_SOURCE_DIR}/../ShapeRecogn/Swig
  ${PROJECT_BINARY_DIR}/doc)

set(medcoupling_LIB_dependancies ${PLATFORM_LIBS} medcouplingremapper medicoco)
if(WIN32 OR CYGWIN)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  list(APPEND medcoupling_LIB_dependancies ${PYTHON_LIBRARIES})
endif()

if(MEDCOUPLING_ENABLE_SHAPERECOGN)
  list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_SHAPE_RECOGN)
  list(APPEND medcoupling_LIB_dependancies shaperecogn)
endif(MEDCOUPLING_ENABLE_SHAPERECOGN)

if(NOT MEDCOUPLING_MICROMED)
  list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_MED_FILE)
  list(APPEND medcoupling_LIB_dependancies medloader)

  if(MEDCOUPLING_ENABLE_RENUMBER)
    list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_RENUMBER)
    list(APPEND medcoupling_LIB_dependancies renumbercpp)
  endif(MEDCOUPLING_ENABLE_RENUMBER)

  if(MEDCOUPLING_ENABLE_PARTITIONER)
    list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_PARTITIONER)
    list(APPEND medcoupling_LIB_dependancies medpartitionercpp)
  endif(MEDCOUPLING_ENABLE_PARTITIONER)
endif(NOT MEDCOUPLING_MICROMED)

if(MEDCOUPLING_USE_MPI)
  include_directories(${MPI_INCLUDE_DIRS})
  add_definitions(${MPI_DEFINITIONS})
  list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_PARALLEL_INTERPOLATOR)
  list(APPEND medcoupling_LIB_dependancies paramedmem)
endif(MEDCOUPLING_USE_MPI)

if(BASICITERATIVESTATISTICS_FOUND)
  list(APPEND SWIG_MODULE_medcoupling_EXTRA_FLAGS -DWITH_ITERATIVE_STATISTICS)
endif(BASICITERATIVESTATISTICS_FOUND)

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(medcoupling python medcoupling.i)
else()
  swig_add_library(
    medcoupling
    LANGUAGE python
    SOURCES medcoupling.i)
endif()

swig_link_libraries(medcoupling ${medcoupling_LIB_dependancies})
swig_check_generation(medcoupling)
if(WIN32 AND NOT MINGW)
  set_target_properties(_medcoupling PROPERTIES DEBUG_OUTPUT_NAME
                                                _medcoupling_d)
  # To increase the size of the .obj file on Windows because
  # MEDLoaderPYTHON_wrap.cxx, generated by SWIG, is too big
  target_compile_options(_medcoupling PRIVATE /bigobj)
endif(WIN32 AND NOT MINGW)

install(TARGETS _medcoupling DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})
install(FILES medcoupling.i medcoupling_pycode
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
salome_install_scripts(
  ${CMAKE_CURRENT_BINARY_DIR}/medcoupling.py ${MEDCOUPLING_INSTALL_PYTHON}
  EXTRA_DPYS "${SWIG_MODULE_medcoupling_REAL_NAME}")

if(MEDCOUPLING_BUILD_PY_TESTS)
  include(tests.set)
endif()

salome_generate_tests_environment(tests_env)

foreach(test ${BASE_TESTS})
  get_filename_component(testname ${test} NAME_WE)
  add_test(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
endforeach()

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/PyWrapping)

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
install(FILES ${BASE_TESTS} DESTINATION ${TEST_INSTALL_DIRECTORY})
install(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
