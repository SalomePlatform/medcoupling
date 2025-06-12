# Copyright (C) 2012-2025  CEA, EDF
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

add_definitions(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set_source_files_properties(ParaMEDMEM.i PROPERTIES CPLUSPLUS ON)
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(ParaMEDMEM.i PROPERTIES SWIG_FLAGS "-py3")
else()
  set_source_files_properties(ParaMEDMEM.i PROPERTIES SWIG_DEFINITIONS
                                                      "-shadow")
endif()
set(SWIG_MODULE_ParaMEDMEM_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
if(MEDCOUPLING_USE_64BIT_IDS)
  string(APPEND SWIG_MODULE_ParaMEDMEM_EXTRA_FLAGS
         ";-DMEDCOUPLING_USE_64BIT_IDS")
endif(MEDCOUPLING_USE_64BIT_IDS)

include_directories(
  ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases)

set(SWIG_MODULE_ParaMEDMEM_EXTRA_DEPS
    ${paramedmem_HEADERS_HXX} ${medloader_HEADERS_HXX}
    ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
    ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX})

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(ParaMEDMEM python ParaMEDMEM.i)
else()
  swig_add_library(
    ParaMEDMEM
    LANGUAGE python
    SOURCES ParaMEDMEM.i)
endif()

swig_link_libraries(ParaMEDMEM paramedmem)
if(WIN32)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  swig_link_libraries(ParaMEDMEM ${PYTHON_LIBRARIES})
endif()
swig_check_generation(ParaMEDMEM)

#
# Tests
#
if(MEDCOUPLING_BUILD_PY_TESTS)
  salome_accumulate_environment(PYTHONPATH NOCHECK
                                ${CMAKE_CURRENT_BINARY_DIR}/../PyWrapping)
  salome_accumulate_environment(MED_RESOURCES_DIR NOCHECK
                                ${CMAKE_SOURCE_DIR}/resources)
  salome_generate_tests_environment(tests_env)

  # -- some tests require 2, 3, 4 or 5 procs -- MPICH does not support
  # --oversubscribe:
  if(NOT ${MPIEXEC_EXECUTABLE} MATCHES "mpich")
    set(_oversub_opt "--oversubscribe")
  endif()

  add_test(NAME PyPara_Basics_Proc2
           COMMAND ${MPIEXEC} -np 2 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_BasicOperation.py)
  set_tests_properties(PyPara_Basics_Proc2 PROPERTIES ENVIRONMENT
                                                      "${tests_env}")

  add_test(NAME PyPara_InterpKernelDEC_Proc4
           COMMAND ${MPIEXEC} -np 4 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_InterpKernelDEC.py)
  set_tests_properties(PyPara_InterpKernelDEC_Proc4 PROPERTIES ENVIRONMENT
                                                               "${tests_env}")
  add_test(NAME PyPara_InterpKernelDEC_Proc5
           COMMAND ${MPIEXEC} -np 5 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_InterpKernelDEC.py)
  set_tests_properties(PyPara_InterpKernelDEC_Proc5 PROPERTIES ENVIRONMENT
                                                               "${tests_env}")

  add_test(NAME PyPara_InterpKernelDEC_easy_Proc5
           COMMAND ${MPIEXEC} -np 5 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_InterpKernelDEC_easy.py)
  set_tests_properties(PyPara_InterpKernelDEC_easy_Proc5
                       PROPERTIES ENVIRONMENT "${tests_env}")

  # ADD_TEST(NAME PyPara_NonCoincidentDEC_Proc5 COMMAND ${MPIEXEC} -np 5
  # ${_oversub_opt} ${PYTHON_EXECUTABLE}
  # ${CMAKE_CURRENT_SOURCE_DIR}/test_NonCoincidentDEC.py)
  # SET_TESTS_PROPERTIES(PyPara_NonCoincidentDEC_Proc5 PROPERTIES ENVIRONMENT
  # "${tests_env}")

  add_test(NAME PyPara_StructuredCoincidentDEC_Proc4
           COMMAND ${MPIEXEC} -np 4 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_StructuredCoincidentDEC.py)
  set_tests_properties(PyPara_StructuredCoincidentDEC_Proc4
                       PROPERTIES ENVIRONMENT "${tests_env}")

  add_test(NAME PyPara_OverlapDEC_Proc4
           COMMAND ${MPIEXEC} -np 4 ${_oversub_opt} ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_OverlapDEC.py)
  set_tests_properties(PyPara_OverlapDEC_Proc4 PROPERTIES ENVIRONMENT
                                                          "${tests_env}")
endif()

set(_tst_scripts test_InterpKernelDEC.py test_NonCoincidentDEC.py
                 test_StructuredCoincidentDEC.py test_OverlapDEC.py)

install(TARGETS _ParaMEDMEM DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})
install(FILES ParaMEDMEM.i ParaMEDMEMCommon.i
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
install(FILES ${_tst_scripts} DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
salome_install_scripts(
  ${CMAKE_CURRENT_BINARY_DIR}/ParaMEDMEM.py ${MEDCOUPLING_INSTALL_PYTHON}
  EXTRA_DPYS "${SWIG_MODULE_ParaMEDMEM_REAL_NAME}")

install(FILES test_InterpKernelDEC.py test_NonCoincidentDEC.py
              test_StructuredCoincidentDEC.py
        DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ParaMEDMEM_Swig)
install(
  FILES test_InterpKernelDEC.py
        test_InterpKernelDEC_easy.py
        test_NonCoincidentDEC.py
        test_OverlapDEC.py
        test_StructuredCoincidentDEC.py
        ParaMEDMEMTestTools.py
        test_BasicOperation.py
  DESTINATION ${TEST_INSTALL_DIRECTORY})
# export MPIEXEC and _oversub_opt to CTestTestfile.cmake of salome test
# mechanism
configure_file(CTestTestfileInstall.cmake.in "CTestTestfileST.cmake" @ONLY)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfileST.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
