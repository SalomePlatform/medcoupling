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

add_definitions(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

set_source_files_properties(MEDRenumber.i PROPERTIES CPLUSPLUS ON)
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(MEDRenumber.i PROPERTIES SWIG_FLAGS "-py3")
else()
  set_source_files_properties(MEDRenumber.i PROPERTIES SWIG_DEFINITIONS
                                                       "-shadow")
endif()
set(SWIG_MODULE_MEDRenumber_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS};-DWITHOUT_AUTOFIELD")
if(MEDCOUPLING_USE_64BIT_IDS)
  string(APPEND SWIG_MODULE_MEDRenumber_EXTRA_FLAGS
         ";-DMEDCOUPLING_USE_64BIT_IDS")
endif(MEDCOUPLING_USE_64BIT_IDS)

set(MEDRenumber_SWIG_DPYS_FILES MEDRenumberCommon.i MEDRenumberImpl.i)

include_directories(
  ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR} # pthread dependency due to python2.7 library
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../RENUMBER
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

if(WIN32)
  set_property(SOURCE MEDRenumber.i PROPERTY COMPILE_DEFINITIONS WIN32)
endif()
if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(MEDRenumber python MEDRenumber.i)
else()
  swig_add_library(
    MEDRenumber
    LANGUAGE python
    SOURCES MEDRenumber.i)
endif()

swig_link_libraries(MEDRenumber ${PLATFORM_LIBS} renumbercpp)
if(WIN32)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  swig_link_libraries(MEDRenumber ${PYTHON_LIBRARIES})
endif()
swig_check_generation(MEDRenumber)

if(WIN32)
  set_target_properties(_MEDRenumber PROPERTIES DEBUG_OUTPUT_NAME
                                                _MEDRenumber_d)
endif(WIN32)
install(TARGETS ${SWIG_MODULE_MEDRenumber_REAL_NAME}
        DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

salome_install_scripts(
  ${CMAKE_CURRENT_BINARY_DIR}/MEDRenumber.py ${MEDCOUPLING_INSTALL_PYTHON}
  EXTRA_DPYS "${SWIG_MODULE_MEDRenumber_REAL_NAME}")

install(FILES MEDRenumber.i MEDRenumberCommon.i MEDRenumberImpl.i
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
install(FILES MEDRenumberTest.py
        DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

if(MEDCOUPLING_BUILD_PY_TESTS)
  salome_accumulate_environment(PYTHONPATH NOCHECK
                                ${CMAKE_CURRENT_BINARY_DIR}/../PyWrapping)
  salome_generate_tests_environment(tests_env)
  add_test(MEDRenumberTest ${PYTHON_EXECUTABLE}
           ${CMAKE_CURRENT_SOURCE_DIR}/MEDRenumberTest.py)
  set_tests_properties(MEDRenumberTest PROPERTIES ENVIRONMENT "${tests_env}")
  if(MEDCOUPLING_ENABLE_PARTITIONER)
    add_test(MEDRenumberUserGuide ${PYTHON_EXECUTABLE}
             ${CMAKE_CURRENT_SOURCE_DIR}/UsersGuideExamplesTest.py)
    set_tests_properties(MEDRenumberUserGuide PROPERTIES ENVIRONMENT
                                                         "${tests_env}")
  endif()
endif()

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/RENUMBER_Swig)
install(FILES MEDRenumberTest.py UsersGuideExamplesTest.py
        DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
