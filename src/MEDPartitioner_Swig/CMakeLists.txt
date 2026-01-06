# Copyright (C) 2012-2026  CEA, EDF
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

include(${SWIG_USE_FILE})

add_definitions(
  ${PYTHON_DEFINITIONS} ${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS}
  ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

set_source_files_properties(MEDPartitioner.i PROPERTIES CPLUSPLUS ON)
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(MEDPartitioner.i PROPERTIES SWIG_FLAGS "-py3")
else()
  set_source_files_properties(MEDPartitioner.i PROPERTIES SWIG_DEFINITIONS
                                                          "-shadow")
endif()
set(SWIG_MODULE_MEDPartitioner_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS};-DWITHOUT_AUTOFIELD")

set(MEDPartitioner_SWIG_DPYS_FILES MEDPartitionerCommon.i)

include_directories(
  ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR} # pthread dependency due to python2.7 library
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDPartitioner
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(MEDPartitioner python MEDPartitioner.i)
else()
  swig_add_library(
    MEDPartitioner
    LANGUAGE python
    SOURCES MEDPartitioner.i)
endif()

swig_link_libraries(MEDPartitioner ${PLATFORM_LIBS} medpartitionercpp)
if(WIN32)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  swig_link_libraries(MEDPartitioner ${PYTHON_LIBRARIES})
endif()
swig_check_generation(MEDPartitioner)

if(WIN32)
  set_target_properties(_MEDPartitioner PROPERTIES DEBUG_OUTPUT_NAME
                                                   _MEDPartitioner_d)
endif(WIN32)
install(TARGETS ${SWIG_MODULE_MEDPartitioner_REAL_NAME}
        DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

set(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDPartitioner.py)
salome_install_scripts("${PYFILES_TO_INSTALL}" ${MEDCOUPLING_INSTALL_PYTHON}
                       EXTRA_DPYS "${SWIG_MODULE_MEDPartitioner_REAL_NAME}")

install(FILES MEDPartitioner.i MEDPartitionerCommon.i
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
install(FILES MEDPartitionerTest.py
        DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

if(MEDCOUPLING_PARTITIONER_METIS)
  if(MEDCOUPLING_BUILD_PY_TESTS)
    # Python test of partitionner require METIS or PARMETIS?
    salome_accumulate_environment(PYTHONPATH NOCHECK
                                  ${CMAKE_CURRENT_BINARY_DIR}/../MEDLoader/Swig)
    salome_generate_tests_environment(tests_env)
    add_test(MEDPartitionerTest ${PYTHON_EXECUTABLE}
             ${CMAKE_CURRENT_SOURCE_DIR}/MEDPartitionerTest.py)
    set_tests_properties(MEDPartitionerTest PROPERTIES ENVIRONMENT
                                                       "${tests_env}")
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader/Swig/MEDLoaderDataForTest.py
      ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  endif()
endif()

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDPartitioner_Swig)
install(FILES MEDPartitionerTest.py DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
