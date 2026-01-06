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
if(MEDCOUPLING_BUILD_PY_TESTS)
  include(tests.set)
endif()

include(${SWIG_USE_FILE})

add_definitions(${PYTHON_DEFINITIONS} ${NUMPY_DEFINITIONS} ${SCIPY_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

set_source_files_properties(MEDCoupling.i PROPERTIES CPLUSPLUS ON)
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(MEDCoupling.i PROPERTIES SWIG_FLAGS "-py3")
else()
  set_source_files_properties(MEDCoupling.i PROPERTIES SWIG_DEFINITIONS
                                                       "-shadow")
endif()
set(SWIG_MODULE_MEDCoupling_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
if(MEDCOUPLING_USE_64BIT_IDS)
  string(APPEND SWIG_MODULE_MEDCoupling_EXTRA_FLAGS
         ";-DMEDCOUPLING_USE_64BIT_IDS")
endif()

set(MEDCoupling_SWIG_DPYS_FILES
    MEDCouplingCommon.i
    MEDCouplingMemArray.i
    DataArrayInt.i
    MEDCouplingFieldDiscretization.i
    MEDCouplingFinalize.i
    MEDCouplingTypemaps.i
    ../ICoCo/Swig/ICoCoMEDField.i)

set(MC_pyTestFiles ${ALL_TESTS} MEDCouplingDataForTest.py)

set(MC_Swig_interf
    MEDCoupling.i
    MEDCouplingCommon.i
    MEDCouplingRemapperCommon.i
    MEDCouplingRefCountObject.i
    MEDCouplingMemArray.i
    DataArrayInt.i
    MEDCouplingFieldDiscretization.i
    MEDCouplingTimeDiscretization.i
    MEDCouplingFinalize.i
    MEDCouplingRemapper.i
    MEDCouplingTypemaps.i
    MEDCouplingDataArrayTypemaps.i
    MEDCouplingDataArrayTraits.hxx
    ../ICoCo/Swig/ICoCoMEDField.i)

include_directories(
  ${NUMPY_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR} # pthread dependency due to python2.7 library
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo/Swig
  ${PROJECT_BINARY_DIR}/doc)

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
set(SWIG_MODULE_MEDCoupling_EXTRA_DEPS
    ${MEDCoupling_SWIG_DPYS_FILES} ${medcoupling_HEADERS_HXX}
    ${medcoupling_HEADERS_TXX} ${interpkernel_HEADERS_HXX}
    ${interpkernel_HEADERS_TXX})

# SWIG must run after the doc if we want to have the docstrings extracted from
# Doxygen into the Python module:
if(MEDCOUPLING_BUILD_DOC)
  list(APPEND SWIG_MODULE_MEDCoupling_EXTRA_FLAGS -DWITH_DOCSTRINGS)
  list(APPEND SWIG_MODULE_MEDCoupling_EXTRA_DEPS
       ${PROJECT_BINARY_DIR}/doc/MEDCoupling_doc.i swig_ready)
endif()

if(WIN32)
  # Under Windows MEDCoupling conflicts wiht medcoupling from PyWrapping
  set(MEDCouling_target_name MEDCouplingCompat)
  set_property(SOURCE MEDCoupling.i PROPERTY SWIG_MODULE_NAME MEDCouplingCompat)
  set_property(SOURCE MEDCoupling.i PROPERTY COMPILE_DEFINITIONS WIN32)
  set(SWIG_MODULE_MEDCouplingCompat_EXTRA_FLAGS
      "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
  if(MEDCOUPLING_USE_64BIT_IDS)
    string(APPEND SWIG_MODULE_MEDCouplingCompat_EXTRA_FLAGS
           ";-DMEDCOUPLING_USE_64BIT_IDS")
  endif(MEDCOUPLING_USE_64BIT_IDS)
else()
  set(MEDCouling_target_name MEDCoupling)
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(${MEDCouling_target_name} python MEDCoupling.i)
else()
  swig_add_library(
    ${MEDCouling_target_name}
    LANGUAGE python
    SOURCES MEDCoupling.i)
endif()

swig_link_libraries(${MEDCouling_target_name} ${PLATFORM_LIBS} medcouplingcpp)
if(WIN32 OR CYGWIN)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  swig_link_libraries(${MEDCouling_target_name} ${PYTHON_LIBRARIES})
endif()

swig_check_generation(${MEDCouling_target_name})
set_source_files_properties(MEDCouplingRemapper.i PROPERTIES CPLUSPLUS ON)
if(WIN32)
  set_property(SOURCE MEDCouplingRemapper.i PROPERTY COMPILE_DEFINITIONS WIN32)
endif()
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  set_source_files_properties(MEDCouplingRemapper.i PROPERTIES SWIG_FLAGS
                                                               "-py3")
else()
  set_source_files_properties(MEDCouplingRemapper.i PROPERTIES SWIG_DEFINITIONS
                                                               "-shadow")
endif()
set(SWIG_MODULE_MEDCouplingRemapper_EXTRA_FLAGS
    "${NUMPY_DEFINITIONS};${SCIPY_DEFINITIONS}")
if(MEDCOUPLING_USE_64BIT_IDS)
  string(APPEND SWIG_MODULE_MEDCouplingRemapper_EXTRA_FLAGS
         ";-DMEDCOUPLING_USE_64BIT_IDS")
endif(MEDCOUPLING_USE_64BIT_IDS)

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
set(SWIG_MODULE_MEDCouplingRemapper_EXTRA_DEPS
    ${MEDCoupling_SWIG_DPYS_FILES} ${medcoupling_HEADERS_HXX}
    ${medcoupling_HEADERS_TXX} ${interpkernel_HEADERS_HXX}
    ${interpkernel_HEADERS_TXX})

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(MEDCouplingRemapper python MEDCouplingRemapper.i)
else()
  swig_add_library(
    MEDCouplingRemapper
    LANGUAGE python
    SOURCES MEDCouplingRemapper.i)
endif()

swig_link_libraries(MEDCouplingRemapper medcouplingremapper)
if(WIN32 OR CYGWIN)
  # PYTHON_LIBRARIES only needed on win32 see
  # https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1
  swig_link_libraries(MEDCouplingRemapper ${PYTHON_LIBRARIES})
endif()

if(WIN32)
  set_target_properties(_MEDCouplingRemapper PROPERTIES DEBUG_OUTPUT_NAME
                                                        _MEDCouplingRemapper_d)
  set_target_properties(
    _${MEDCouling_target_name} PROPERTIES DEBUG_OUTPUT_NAME
                                          _${MEDCouling_target_name}_d)
endif(WIN32)
install(TARGETS ${SWIG_MODULE_${MEDCouling_target_name}_REAL_NAME}
                ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}
        DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

set(PYFILES_TO_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingRemapper.py)
if(WIN32)
  set(PYFILES_TO_INSTALL ${PYFILES_TO_INSTALL}
                         ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingCompat.py)
else()
  set(PYFILES_TO_INSTALL ${PYFILES_TO_INSTALL}
                         ${CMAKE_CURRENT_BINARY_DIR}/MEDCoupling.py)
endif()
install(FILES vtk2medcoupling.py geom2medcoupling.py
        DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

salome_install_scripts(
  "${PYFILES_TO_INSTALL}"
  ${MEDCOUPLING_INSTALL_PYTHON}
  EXTRA_DPYS
  "${SWIG_MODULE_${MEDCouling_target_name}_REAL_NAME};${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}"
)

install(FILES ${MC_Swig_interf} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
install(FILES ${MC_pyTestFiles}
        DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})

salome_accumulate_environment(PYTHONPATH NOCHECK
                              ${CMAKE_CURRENT_BINARY_DIR}/../PyWrapping)
salome_generate_tests_environment(tests_env)

foreach(test ${BASE_TESTS})
  get_filename_component(testname ${test} NAME_WE)
  add_test(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
endforeach()

if(NUMPY_FOUND)
  foreach(test ${NUMPY_TESTS})
    get_filename_component(testname ${test} NAME_WE)
    add_test(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE}
                                      ${CMAKE_CURRENT_SOURCE_DIR}/${test})
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
  endforeach()
endif(NUMPY_FOUND)

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDCoupling_Swig)
install(FILES ${MC_pyTestFiles} DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
install(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
