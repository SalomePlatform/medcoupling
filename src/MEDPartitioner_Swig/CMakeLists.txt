# Copyright (C) 2012-2024  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

INCLUDE(${SWIG_USE_FILE})

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

IF (NOT DEFINED MSVC)
  ADD_DEFINITIONS(-Wsign-compare -Wconversion)
ENDIF()

SET_SOURCE_FILES_PROPERTIES(MEDPartitioner.i PROPERTIES CPLUSPLUS ON)
IF ("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  SET_SOURCE_FILES_PROPERTIES(MEDPartitioner.i PROPERTIES SWIG_FLAGS "-py3")
ELSE()
  SET_SOURCE_FILES_PROPERTIES(MEDPartitioner.i PROPERTIES SWIG_DEFINITIONS "-shadow")
ENDIF()
SET(SWIG_MODULE_MEDPartitioner_EXTRA_FLAGS "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS};-DWITHOUT_AUTOFIELD")

SET (MEDPartitioner_SWIG_DPYS_FILES
    MEDPartitionerCommon.i)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR} # pthread dependency due to python2.7 library
  ${NUMPY_INCLUDE_DIR}
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
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

IF(${CMAKE_VERSION} VERSION_LESS "3.8.0")     
  SWIG_ADD_MODULE(MEDPartitioner python MEDPartitioner.i)
ELSE()
  SWIG_ADD_LIBRARY(MEDPartitioner LANGUAGE python SOURCES MEDPartitioner.i)
ENDIF()

SWIG_LINK_LIBRARIES(MEDPartitioner ${PYTHON_LIBRARIES} ${PLATFORM_LIBS} medpartitionercpp)
SWIG_CHECK_GENERATION(MEDPartitioner)

IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDPartitioner PROPERTIES DEBUG_OUTPUT_NAME _MEDPartitioner_d)
ENDIF(WIN32)
INSTALL(TARGETS ${SWIG_MODULE_MEDPartitioner_REAL_NAME} DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDPartitioner.py)
SALOME_INSTALL_SCRIPTS("${PYFILES_TO_INSTALL}" ${MEDCOUPLING_INSTALL_PYTHON}  EXTRA_DPYS "${SWIG_MODULE_MEDPartitioner_REAL_NAME}")

INSTALL(FILES MEDPartitioner.i MEDPartitionerCommon.i DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
INSTALL(FILES MEDPartitionerTest.py DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

IF(MEDCOUPLING_PARTITIONER_METIS)
    IF(MEDCOUPLING_BUILD_PY_TESTS)
       # Python test of partitionner require METIS or PARMETIS?
       SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK ${CMAKE_CURRENT_BINARY_DIR}/../MEDLoader/Swig)
       SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
       ADD_TEST(MEDPartitionerTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDPartitionerTest.py)
       SET_TESTS_PROPERTIES(MEDPartitionerTest PROPERTIES ENVIRONMENT "${tests_env}")
       CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader/Swig/MEDLoaderDataForTest.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
    ENDIF()
ENDIF()

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDPartitioner_Swig)
INSTALL(FILES MEDPartitionerTest.py DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
