# Copyright (C) 2012-2024  CEA, EDF
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

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

set_source_files_properties(MEDCoupling.i PROPERTIES CPLUSPLUS ON)

# set(MC_pyTestFiles ${ALL_TESTS} MEDCouplingDataForTest.py)

# _ABR_ Ensure dependency mechanism on all SWIG files and headers
# set(SWIG_MODULE_MEDCoupling_EXTRA_DEPS
#     ${MEDCoupling_SWIG_DPYS_FILES} ${medcoupling_HEADERS_HXX}
#     ${medcoupling_HEADERS_TXX} ${interpkernel_HEADERS_HXX}
#     ${interpkernel_HEADERS_TXX})

# SWIG must run after the doc if we want to have the docstrings extracted from
# Doxygen into the Python module:
# if(MEDCOUPLING_BUILD_DOC)
#   list(APPEND SWIG_MODULE_MEDCoupling_EXTRA_FLAGS -DWITH_DOCSTRINGS)
#   list(APPEND SWIG_MODULE_MEDCoupling_EXTRA_DEPS
#        ${PROJECT_BINARY_DIR}/doc/MEDCoupling_doc.i swig_ready)
# endif()

set_property(SOURCE MEDCoupling.i PROPERTY SWIG_MODULE_NAME pymedcoupling)

if(WIN32)
  # Under Windows MEDCoupling conflicts wiht medcoupling from PyWrapping
  set_property(SOURCE MEDCoupling.i PROPERTY COMPILE_DEFINITIONS WIN32)
endif()

swig_add_library(
  pymedcoupling
  LANGUAGE python
  OUTPUT_DIR .
  SOURCES MEDCoupling.i)

target_include_directories(
  ${SWIG_MODULE_pymedcoupling_REAL_NAME}
  PRIVATE
  "\$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "\$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../icoco>")

set_property(TARGET ${SWIG_MODULE_pymedcoupling_REAL_NAME} PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
set_property(TARGET ${SWIG_MODULE_pymedcoupling_REAL_NAME} PROPERTY INSTALL_RPATH ${MEDCOUPLING_INSTALL_PYTHON_LIBS_DEPS_RPATH})

target_link_libraries(${SWIG_MODULE_pymedcoupling_REAL_NAME}
  PUBLIC medcouplingcpp
  PRIVATE Python::Module Python::NumPy)

###################
# Remapper target #
###################

set_source_files_properties(MEDCouplingRemapper.i PROPERTIES CPLUSPLUS ON)
if(WIN32)
  set_property(SOURCE MEDCouplingRemapper.i PROPERTY COMPILE_DEFINITIONS WIN32)
endif()

swig_add_library(
  MEDCouplingRemapper
  LANGUAGE python
  OUTPUT_DIR .
  SOURCES MEDCouplingRemapper.i)

target_include_directories(
  ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}
  PRIVATE
  "\$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "\$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../icoco>")

set_property(TARGET ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME} PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
set_property(TARGET ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME} PROPERTY INSTALL_RPATH ${MEDCOUPLING_INSTALL_PYTHON_LIBS_DEPS_RPATH})

target_link_libraries(${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}
  PUBLIC medcouplingremapper
  PRIVATE Python::Module Python::NumPy)

if(WIN32)
  set_target_properties(${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}
    PROPERTIES DEBUG_OUTPUT_NAME _MEDCouplingRemapper_d)
  set_target_properties(
    ${SWIG_MODULE_pymedcoupling_REAL_NAME}
    PROPERTIES DEBUG_OUTPUT_NAME _${MEDCouling_target_name}_d)
endif()

################
# Installation #
################

install(TARGETS ${SWIG_MODULE_pymedcoupling_REAL_NAME}
  ${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}
  # RUNTIME_DEPENDENCIES DIRECTORIES ${MEDCOUPLING_INSTALL_LIBS}
  DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/MEDCouplingRemapper.py
  DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/pymedcoupling.py
  DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

install(FILES vtk2medcoupling.py geom2medcoupling.py
  DESTINATION ${MEDCOUPLING_INSTALL_PYTHON})

# salome_install_scripts(
#   "${PYFILES_TO_INSTALL}"
#   ${MEDCOUPLING_INSTALL_PYTHON}
#   EXTRA_DPYS
#   "${SWIG_MODULE_${MEDCouling_target_name}_REAL_NAME};${SWIG_MODULE_MEDCouplingRemapper_REAL_NAME}"
# )
# 
# install(FILES ${MC_Swig_interf} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
# install(FILES ${MC_pyTestFiles}
#         DESTINATION ${MEDCOUPLING_INSTALL_SCRIPT_PYTHON})
# 
# salome_accumulate_environment(PYTHONPATH NOCHECK
#                               ${CMAKE_CURRENT_BINARY_DIR}/../PyWrapping)
# salome_generate_tests_environment(tests_env)

set(BASE_TESTS
  MEDCouplingBasicsTest1.py
  MEDCouplingBasicsTest2.py
  MEDCouplingBasicsTest3.py
  MEDCouplingBasicsTest4.py
  MEDCouplingBasicsTest5.py
  MEDCouplingBasicsTest6.py
  MEDCouplingBasicsTest7.py
  MEDCouplingIntersectTest.py
  MEDCouplingExamplesTest.py
  MEDCouplingRemapperTest.py
  UsersGuideExamplesTest.py
)

# if numpy is used
set(NUMPY_TESTS
  MEDCouplingNumPyTest.py
  MEDCouplingPickleTest.py
  UsersGuideExamplesTest_numpy.py
)

set(ALL_TESTS ${BASE_TESTS} ${NUMPY_TESTS})

foreach(test ${BASE_TESTS})
  get_filename_component(testname ${test} NAME_WE)
  add_test(NAME ${testname} COMMAND ${Python_EXECUTABLE}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
endforeach()

if(NUMPY_FOUND)
  foreach(test ${NUMPY_TESTS})
    get_filename_component(testname ${test} NAME_WE)
    add_test(NAME ${testname} COMMAND ${Python_EXECUTABLE}
                                      ${CMAKE_CURRENT_SOURCE_DIR}/${test})
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
  endforeach()
endif(NUMPY_FOUND)

# # Application tests
# 
# set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDCoupling_Swig)
# install(FILES ${MC_pyTestFiles} DESTINATION ${TEST_INSTALL_DIRECTORY})
# 
# install(
#   FILES CTestTestfileInstall.cmake
#   DESTINATION ${TEST_INSTALL_DIRECTORY}
#   RENAME CTestTestfile.cmake)
# install(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
