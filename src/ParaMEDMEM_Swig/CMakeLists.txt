# Copyright (C) 2012-2021  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

IF (NOT DEFINED MSVC)
  ADD_DEFINITIONS(-Wsign-compare -Wconversion)
ENDIF()

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
IF ("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-py3")
ELSE()
  SET_SOURCE_FILES_PROPERTIES(ParaMEDMEM.i PROPERTIES SWIG_DEFINITIONS "-shadow")
ENDIF()
SET(SWIG_MODULE_ParaMEDMEM_EXTRA_FLAGS "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
IF(MEDCOUPLING_USE_64BIT_IDS)
  STRING(APPEND SWIG_MODULE_ParaMEDMEM_EXTRA_FLAGS ";-DMEDCOUPLING_USE_64BIT_IDS")
ENDIF(MEDCOUPLING_USE_64BIT_IDS)

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIR}
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  )

SET (SWIG_MODULE_ParaMEDMEM_EXTRA_DEPS
    ${paramedmem_HEADERS_HXX}
    ${medloader_HEADERS_HXX}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})

IF(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  SWIG_ADD_MODULE(ParaMEDMEM python ParaMEDMEM.i)
ELSE()
  SWIG_ADD_LIBRARY(ParaMEDMEM LANGUAGE python SOURCES ParaMEDMEM.i)
ENDIF()

SWIG_LINK_LIBRARIES(ParaMEDMEM ${PYTHON_LIBRARIES} paramedmem)
SWIG_CHECK_GENERATION(ParaMEDMEM)

#
# Tests
#
IF(MEDCOUPLING_BUILD_PY_TESTS)
    SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK ${CMAKE_CURRENT_BINARY_DIR}/../PyWrapping)
    SALOME_ACCUMULATE_ENVIRONMENT(MED_RESOURCES_DIR NOCHECK ${CMAKE_SOURCE_DIR}/resources)
    SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
    
    #  -- some tests require 2, 3, 4 or 5 procs --
    # MPICH does not support --oversubscribe:
    IF(NOT ${MPIEXEC_EXECUTABLE} MATCHES "mpich")
        SET(_oversub_opt "--oversubscribe")
    ENDIF()
    
    ADD_TEST(NAME PyPara_Basics_Proc2
             COMMAND ${MPIEXEC} -np 2  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_BasicOperation.py)
    SET_TESTS_PROPERTIES(PyPara_Basics_Proc2 PROPERTIES ENVIRONMENT "${tests_env}")
    
    ADD_TEST(NAME PyPara_InterpKernelDEC_Proc4
             COMMAND ${MPIEXEC} -np 4  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_InterpKernelDEC.py)
    SET_TESTS_PROPERTIES(PyPara_InterpKernelDEC_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")
    ADD_TEST(NAME PyPara_InterpKernelDEC_Proc5
             COMMAND ${MPIEXEC} -np 5  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_InterpKernelDEC.py)
    SET_TESTS_PROPERTIES(PyPara_InterpKernelDEC_Proc5 PROPERTIES ENVIRONMENT "${tests_env}")
    
    #ADD_TEST(NAME PyPara_NonCoincidentDEC_Proc5
    #         COMMAND ${MPIEXEC} -np 5  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_NonCoincidentDEC.py)
    #SET_TESTS_PROPERTIES(PyPara_NonCoincidentDEC_Proc5 PROPERTIES ENVIRONMENT "${tests_env}")
    
    ADD_TEST(NAME PyPara_StructuredCoincidentDEC_Proc4
             COMMAND ${MPIEXEC} -np 4  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_StructuredCoincidentDEC.py)
    SET_TESTS_PROPERTIES(PyPara_StructuredCoincidentDEC_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")
    
    ADD_TEST(NAME PyPara_OverlapDEC_Proc4
             COMMAND ${MPIEXEC} -np 4  ${_oversub_opt} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_OverlapDEC.py)
    SET_TESTS_PROPERTIES(PyPara_OverlapDEC_Proc4 PROPERTIES ENVIRONMENT "${tests_env}")
ENDIF()

SET(_tst_scripts
  test_InterpKernelDEC.py
  test_NonCoincidentDEC.py
  test_StructuredCoincidentDEC.py
  test_OverlapDEC.py
)

INSTALL(TARGETS _ParaMEDMEM DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})
INSTALL(FILES ParaMEDMEM.i ParaMEDMEMCommon.i DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
INSTALL(FILES ${_tst_scripts} DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
SALOME_INSTALL_SCRIPTS(${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py ${MEDCOUPLING_INSTALL_PYTHON} EXTRA_DPYS "${SWIG_MODULE_ParaMEDMEM_REAL_NAME}")

INSTALL(FILES test_InterpKernelDEC.py test_NonCoincidentDEC.py test_StructuredCoincidentDEC.py DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ParaMEDMEM_Swig)
install(FILES test_InterpKernelDEC.py test_NonCoincidentDEC.py test_OverlapDEC.py test_StructuredCoincidentDEC.py ParaMEDMEMTestTools.py test_BasicOperation.py DESTINATION ${TEST_INSTALL_DIRECTORY})
# export MPIEXEC and _oversub_opt to CTestTestfile.cmake of salome test mechanism
configure_file(CTestTestfileInstall.cmake.in "CTestTestfileST.cmake" @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfileST.cmake DESTINATION ${TEST_INSTALL_DIRECTORY} RENAME CTestTestfile.cmake)
