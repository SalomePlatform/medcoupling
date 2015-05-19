# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${PYTHON_DEFINITIONS})

SET_SOURCE_FILES_PROPERTIES(MEDPartitioner.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDPartitioner.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDPartitioner_EXTRA_FLAGS "${SWIG_EXTRA_FLAGS_FOR_NUMPYANDSCIPY};-DWITHOUT_AUTOFIELD")

SET (MEDPartitioner_SWIG_DPYS_FILES
    MEDPartitionerCommon.i)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR} # pthread dependancy due to python2.7 library
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
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

SWIG_ADD_MODULE(MEDPartitioner python MEDPartitioner.i)
SWIG_LINK_LIBRARIES(MEDPartitioner ${PYTHON_LIBRARIES} ${PLATFORM_LIBS} medpartitionercpp)

IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDPartitioner PROPERTIES DEBUG_OUTPUT_NAME _MEDPartitioner_d)
ENDIF(WIN32)
INSTALL(TARGETS ${SWIG_MODULE_MEDPartitioner_REAL_NAME} DESTINATION ${SALOME_INSTALL_PYTHON})

SET(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDPartitioner.py)
INSTALL_AND_COMPILE_PYTHON_FILE("${PYFILES_TO_INSTALL}" ${SALOME_INSTALL_SCRIPT_PYTHON})

INSTALL(FILES MEDPartitioner.i MEDPartitionerCommon.i DESTINATION ${SALOME_INSTALL_HEADERS})
INSTALL(FILES MEDPartitionerTest.py DESTINATION ${SALOME_INSTALL_SCRIPT_PYTHON})

SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
ADD_TEST(MEDPartitionerTest ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MEDPartitionerTest.py)
SET_TESTS_PROPERTIES(MEDPartitionerTest PROPERTIES ENVIRONMENT "${tests_env}")

# Application tests

SET(TEST_INSTALL_DIRECTORY ${SALOME_INSTALL_SCRIPT_SCRIPTS}/test/MEDCoupling/MEDPartitioner_Swig)
INSTALL(FILES MEDPartitionerTest.py DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)


